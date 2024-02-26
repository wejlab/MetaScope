globalVariables("count")

obtain_reads <- function(input_file, input_type, aligner, blast_fastas = FALSE, quiet) {
  if (blast_fastas) {
    to_pull <- c("qname", "rname", "cigar", "qwidth", "pos", "seq")
  } else {
    to_pull <- c("qname", "rname", "cigar", "qwidth", "pos")
  }
  if (identical(input_type, "bam")) {
    if (!quiet) message("Reading .bam file: ", input_file)
    if (identical(aligner, "bowtie2")) {
      params <- Rsamtools::ScanBamParam(what = to_pull, tag = c("AS"))
    } else if (identical(aligner, "subread")) {
      params <- Rsamtools::ScanBamParam(what = to_pull, tag = c("NM"))
    } else if (identical(aligner, "other")) {
      params <- Rsamtools::ScanBamParam(what = to_pull)
    }
    reads <- Rsamtools::scanBam(input_file, param = params)
  } else if (identical(input_type, "csv.gz")) {
    if (!quiet) message("Reading .csv.gz file: ", input_file)
    reads <- data.table::fread(input_file, sep = ",", header = FALSE) %>%
      magrittr::set_colnames(c(to_pull, "tag")) %>% as.list() %>% list()
    if (identical(aligner, "bowtie2")) {
      reads[[1]]$tag <- list("AS" = reads[[1]]$tag)
    } else if (identical(aligner, "subread")) {
      reads[[1]]$tag <- list("NM" = reads[[1]]$tag)
    }
  }
  return(reads)
}

identify_rnames <- function(reads, unmapped = NULL) {
  reads_in <- reads[[1]]$rname
  if(!is.null(unmapped)) reads_in <- reads[[1]]$rname[!unmapped]
  # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole
  prefixes <- c("AC", "NC", "NG", "NT", "NW", "NZ") %>%
    paste0("_") %>%
    paste(collapse = "|")
  this_read <- stringr::str_split_i(reads_in, prefixes, -1) |>
    stringr::str_split_i("[| ]", 1)
  read_prefix <- stringr::str_extract(reads_in, prefixes)
  mapped_rname <- paste0(read_prefix, this_read)
  return(mapped_rname)
}

find_accessions <- function(accessions, NCBI_key, quiet) {
  # Convert accessions to taxids and get genome names
  if (!quiet) message("Obtaining taxonomy and genome names")
  # If URI length is greater than 2500 characters then split accession list
  URI_length <- nchar(paste(accessions, collapse = "+"))
  if (URI_length > 2500) {
    chunks <- split(accessions, ceiling(seq_along(accessions) / 100))
    tax_id_all <- c()
    if (!quiet) message("Accession list broken into ", length(chunks),
                        " chunks")
    for (i in seq_along(chunks)) {
      success <- FALSE
      attempt <- 0
      # Attempt to get taxid up to three times for each chunk
      while (!success) {
        try({
          attempt <- attempt + 1
          if (attempt > 1 && !quiet) message(
            "Attempt #", attempt, " Chunk #", i)
          tax_id_chunk <- taxize::genbank2uid(id = chunks[[i]],
                                              key = NCBI_key)
          Sys.sleep(1)
          tax_id_all <- c(tax_id_all, tax_id_chunk)
          success <- TRUE
        })
      }
    }
  } else tax_id_all <- taxize::genbank2uid(id = accessions, key = NCBI_key)
  return(tax_id_all)
}

missing_acc_finder <- function(accessions, NCBI_key = NULL) {
  na_ind <- which(is.na(accessions))
  if(length(na_ind) != 0) {
    result <- accessions[na_ind] %>%
      find_accessions(quiet = TRUE, NCBI_key = NCBI_key) %>%
      plyr::aaply(1, function(x) x[1]) %>% unname()
    accessions[na_ind] <- result
  }
  return(accessions)
}

get_alignscore <- function(aligner, cigar_strings, count_matches, scores,
                           qwidths) {
  #Subread alignment scores: CIGAR string matches - edit score
  if (identical(aligner, "subread")) {
    num_match <- unlist(vapply(cigar_strings, count_matches,
                               USE.NAMES = FALSE, double(1)))
    alignment_scores <- num_match - scores
    scaling_factor <- 100.0 / max(alignment_scores)
    relative_alignment_scores <- alignment_scores - min(alignment_scores)
    exp_alignment_scores <- exp(relative_alignment_scores * scaling_factor)
  } else if (identical(aligner, "bowtie2")) {
    # Bowtie2 alignment scores: AS value + read length (qwidths)
    alignment_scores <- scores + qwidths
    scaling_factor <- 100.0 / max(alignment_scores)
    relative_alignment_scores <- alignment_scores - min(alignment_scores)
    exp_alignment_scores <- exp(relative_alignment_scores * scaling_factor)
  } else if (identical(aligner, "other")) {
    # Other alignment scores: No assumptions
    exp_alignment_scores <- 1
  }
  return(exp_alignment_scores)
}

get_assignments <- function(combined, convEM, maxitsEM, unique_taxids,
                            unique_genome_names, update_bam = TRUE, input_file, quiet) {
  combined$index <- seq.int(1, nrow(combined))
  input_distinct <- dplyr::distinct(combined, .data$qname, .data$rname,
                                    .keep_all = TRUE)
  qname_inds_2 <- input_distinct$qname
  rname_tax_inds_2 <- input_distinct$rname
  scores_2 <- input_distinct$scores
  non_unique_read_ind <- unique(combined[[1]][(
    duplicated(input_distinct[, 1]) | duplicated(input_distinct[, 1],
                                                 fromLast = TRUE))])
  # 1 if read is multimapping, 0 if read is unique
  y_ind_2 <- as.numeric(unique(input_distinct[[1]]) %in% non_unique_read_ind)
  gammas <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2, x = scores_2)
  pi_old <- rep(1 / nrow(gammas), ncol(gammas))
  pi_new <- theta_new <- Matrix::colMeans(gammas)
  conv <- max(abs(pi_new - pi_old) / pi_old)
  it <- 0
  if (!quiet) message("Starting EM iterations")
  while (conv > convEM && it < maxitsEM) {
    # Expectation Step: Estimate expected value for each read to ea genome
    pi_mat <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2,
                                   x = pi_new[rname_tax_inds_2])
    theta_mat <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2,
                                      x = theta_new[rname_tax_inds_2])
    weighted_gamma <- gammas * pi_mat * theta_mat
    weighted_gamma_sums <- Matrix::rowSums(weighted_gamma)
    gammas_new <- weighted_gamma / weighted_gamma_sums
    # Maximization step: proportion of reads to each genome
    pi_new <- Matrix::colMeans(gammas_new)
    theta_new_num <- (Matrix::colSums(y_ind_2 * gammas_new) + 1)
    theta_new <- theta_new_num / (nrow(gammas_new) + 1)
    # Check convergence
    it <- it + 1
    conv <- max(abs(pi_new - pi_old) / pi_old, na.rm = TRUE)
    pi_old <- pi_new
    if (!quiet) message(c(it, conv))
  }
  if (!quiet) message("\tDONE! Converged in ", it, " iterations.")

  hit_which <- qlcMatrix::rowMax(gammas_new, which = TRUE)$which
  hit <- mapply(function(q, r) hit_which[q,r], combined$qname, combined$rname)
  combined$hit <- hit


  combined_single <- combined %>% group_by(qname, rname) %>%
    mutate(best_hit = (hit & max(scores) == scores)) %>%
    ungroup() %>% group_by(qname) %>%
    mutate(single_hit = row_number() == min(row_number()[best_hit]))

  combined_distinct <- dplyr::distinct(combined, .data$qname, .data$rname,
                                        .keep_all = TRUE)
  combined_distinct <- combined_distinct[combined_distinct$hit == TRUE,]
  best_hit <- Matrix::colSums(hit_which)
  names(best_hit) <- seq_along(best_hit)
  best_hit <- best_hit[best_hit != 0]
  hits_ind <- as.numeric(names(best_hit))
  final_taxids <- unique_taxids[hits_ind]
  final_genomes <- unique_genome_names[hits_ind]
  proportion <- best_hit / sum(best_hit)
  gammasums <- Matrix::colSums(gammas_new)
  readsEM <- round(gammasums[hits_ind], 1)
  propEM <- gammasums[hits_ind] / sum(gammas_new)
  results_tibble <- dplyr::tibble(TaxonomyID = final_taxids, Genome = final_genomes,
                           read_count = best_hit, Proportion = proportion,
                           readsEM = readsEM, EMProportion = propEM,
                           hits_ind = hits_ind) %>%
    dplyr::arrange(dplyr::desc(.data$read_count))
  if (!quiet) message("Found reads for ", nrow(results_tibble), " genomes")

  return(list(results_tibble, combined_distinct, combined_single))
}

#' Count the number of base lengths in a CIGAR string for a given operation
#'
#' The 'CIGAR' (Compact Idiosyncratic Gapped Alignment Report) string is how the
#' SAM/BAM format represents spliced alignments. This function will accept a
#' CIGAR string for a single read and a single character indicating the
#' operation to be parsed in the string. An operation is a type of column that
#' appears in the alignment, e.g. a match or gap. The integer following the
#' operator specifies a number of consecutive operations. The
#' \code{count_matches()} function will identify all occurrences of the operator
#' in the string input, add them, and return an integer number representing the
#' total number of operations for the read that was summarized by the input
#' CIGAR string.
#'
#' This function is best used on a vector of CIGAR strings using an apply
#' function (see examples).
#'
#' @param x Character. A CIGAR string for a read to be parsed. Examples of
#' possible operators include "M", "D", "I", "S", "H", "=", "P", and "X".
#' @param char A single letter representing the operation to total for the
#' given string.
#'
#' @return an integer number representing the total number of alignment
#' operations for the read that was summarized by the input CIGAR string.
#'
#' @export
#'
#' @examples
#' # A single cigar string: 3M + 3M + 5M
#' cigar1 <- "3M1I3M1D5M"
#' count_matches(cigar1, char = "M")
#'
#' # Parse with operator "P": 2P
#' cigar2 <- "4M1I2P9M"
#' count_matches(cigar2, char = "P")
#'
#' # Apply to multiple strings: 1I + 1I + 5I
#' cigar3 <- c("3M1I3M1D5M", "4M1I1P9M", "76M13M5I")
#' vapply(cigar3, count_matches, char = "I",
#'        FUN.VALUE = numeric(1))
#'

count_matches <- function(x, char = "M") {
  if (length(char) != 1) {
    stop("Please provide a single character ",
         "operator with which to parse.")
  } else if (length(x) != 1) {
    stop("Please provide a single CIGAR string to be parsed.")
  }
  pattern <- paste0("\\d+", char)
  ind <- gregexpr(pattern, x)[[1]]
  start <- as.numeric(ind)
  end <- start + attr(ind, "match.length") - 2
  out <- cbind(start, end) %>% apply(
    1, function(y) substr(x, start = y[1], stop = y[2])) %>%
    as.numeric() %>% sum()
  return(data.table::fifelse(is.na(out[1]), yes = 0, no = out[1]))
}

#' Helper Function for MetaScope ID
#'
#' Used to create plots of genome coverage for any number of accession numbers
#'
#' @param which_taxid Which taxid to plot
#' @param which_genome Which genome to plot
#' @param accessions List of accessions from \code{metascope_id()}
#' @param taxids List of accessions from \code{metascope_id()}
#' @param reads List of reads from input file
#' @param out_base The basename of the input file
#' @param out_dir The path to the input file
#'
#' @return A plot of the read coverage for a given genome

locations <- function(which_taxid, which_genome,
                      accessions, taxids, reads, out_base, out_dir) {
  plots_save <- file.path(out_dir, paste(out_base, "cov_plots",
                                         sep = "_"))
  # map back to accessions
  choose_acc <- paste(accessions[which(as.numeric(taxids) %in% which_taxid)])
  # map back to BAM
  map2bam_acc <- which(reads[[1]]$rname %in% choose_acc)
  # Split genome name to make digestible
  this_genome <- strsplit(which_genome, " ")[[1]][c(1, 2)]
  use_name <- paste(this_genome, collapse = " ") %>% stringr::str_replace(",", "")
  coverage <- round(mean(seq_len(338099) %in% unique(
    reads[[1]]$pos[map2bam_acc])), 3)
  # Plotting
  dfplot <- dplyr::tibble(x = reads[[1]]$pos[map2bam_acc])
  ggplot2::ggplot(dfplot, ggplot2::aes(.data$x, fill = ggplot2::after_stat(count))) +
    ggplot2::geom_histogram(bins = 50) +
    ggplot2::scale_fill_gradient(low = 'red', high = 'yellow') +
    ggplot2::theme_classic() +
    ggplot2::labs(title = bquote("Positions of reads mapped to"~italic(.(use_name))),
                  xlab = "Aligned position across genome (leftmost read position)",
                  ylab = "Read Count",
                  caption = paste0("Accession Number: ", choose_acc))

  ggplot2::ggsave(paste0(plots_save, "/",
                         stringr::str_replace(use_name, " ", "_"), ".png"),
                  device = "png")
}

#' Identify which genomes are represented in a processed sample
#'
#' This function will read in a .bam or .csv.gz file, annotate the taxonomy and
#' genome names, reduce the mapping ambiguity using a mixture model, and output
#' a .csv file with the results. Currently, it assumes that the genome
#' library/.bam files use NCBI accession names for reference names (rnames in
#' .bam file).
#'
#' @param input_file The .bam or .csv.gz file of sample reads to be identified.
#' @param input_type Extension of file input. Should be either "bam" or
#'   "csv.gz". Default is "csv.gz".
#' @param aligner The aligner which was used to create the bam file. Default is
#'   "bowtie2" but can also be set to "subread" or "other".
#' @param db Currently accepts one of \code{c("ncbi", "silva", "other")}
#' Default is \code{"ncbi"}, appropriate for samples aligned against indices
#' compiled from NCBI whole genome databases. Alternatively, usage of an
#' alternate database (like Greengenes2) should be specified with
#' \code{"other"}.
#' @param db_feature_table If \code{db = "other"}, a data.frame must be supplied
#' with two columns, "Feature ID" matching the names of the alignment indices,
#' and a second \code{character} column supplying the taxon identifying information.
#' @param NCBI_key (character) NCBI Entrez API key. optional. See
#'   taxize::use_entrez(). Due to the high number of requests made to NCBI, this
#'   function will be less prone to errors if you obtain an NCBI key. You may
#'   enter the string as an input or set it as ENTREZ_KEY in .Renviron.
#' @param out_dir The directory to which the .csv output file will be output.
#'   Defaults to \code{dirname(input_file)}.
#' @param convEM The convergence parameter of the EM algorithm. Default set at
#'   \code{1/10000}.
#' @param maxitsEM The maximum number of EM iterations, regardless of whether
#'   the convEM is below the threshhold. Default set at \code{50}. If set at
#'   \code{0}, the algorithm skips the EM step and summarizes the .bam file 'as
#'   is'
#' @param blast_fastas Logical whether or not to output fasta files for MetaBlast
#' @param num_genomes Number of genomes to output fasta files for MetaBlast
#' @param num_reads Number of reads per genome per fasta file for MetaBlast
#' @param num_species_plot The number of genome coverage plots to be saved.
#'   Default is \code{NULL}, which saves coverage plots for the ten most highly
#'   abundant species.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#'
#' @return This function returns a .csv file with annotated read counts to
#'   genomes with mapped reads. The function itself returns the output .csv file
#'   name.
#'
#' @export
#'
#' @examples
#' #### Align reads to reference library and then apply metascope_id()
#' ## Assuming filtered bam files already exist
#'
#' ## Create temporary directory
#' file_temp <- tempfile()
#' dir.create(file_temp)
#'
#' #### Subread aligned bam file
#'
#' ## Create object with path to filtered subread csv.gz file
#' filt_file <- "subread_target.filtered.csv.gz"
#' bamPath <- system.file("extdata", filt_file, package = "MetaScope")
#' file.copy(bamPath, file_temp)
#'
#' ## Run metascope id with the aligner option set to subread
#' metascope_id(input_file = file.path(file_temp, filt_file),
#'              aligner = "subread", num_species_plot = 0,
#'              input_type = "csv.gz")
#'
#' #### Bowtie 2 aligned .csv.gz file
#'
#' ## Create object with path to filtered bowtie2 bam file
#' bowtie_file <- "bowtie_target.filtered.csv.gz"
#' bamPath <- system.file("extdata", bowtie_file, package = "MetaScope")
#' file.copy(bamPath, file_temp)
#'
#' ## Run metascope id with the aligner option set to bowtie2
#' metascope_id(file.path(file_temp, bowtie_file), aligner = "bowtie2",
#'              num_species_plot = 0, input_type = "csv.gz")
#'
#' ## Remove temporary directory
#' unlink(file_temp, recursive = TRUE)
#'

metascope_id <- function(input_file, input_type = "csv.gz",
                         aligner = "bowtie2",
                         db = "ncbi",
                         db_feature_table = NULL,
                         NCBI_key = NULL,
                         out_dir = dirname(input_file),
                         convEM = 1 / 10000, maxitsEM = 25,
                         blast_fastas = TRUE, num_genomes = 100,
                         num_reads = 50, update_bam = TRUE,
                         num_species_plot = NULL,
                         quiet = TRUE)  {
  out_base <- input_file %>% base::basename() %>% strsplit(split = "\\.") %>%
    magrittr::extract2(1) %>% magrittr::extract(1)
  out_file <- file.path(out_dir, paste0(out_base, ".metascope_id.csv"))
  # Check to make sure valid aligner is specified
  if (!(aligner %in% c("bowtie2", "subread", "other"))) {
    stop("Please make sure aligner is set to either 'bowtie2', 'subread',",
         " or 'other'")
  }
  if (db == "other" && is.null(db_feature_table)) {
    stop("Please supply a data.frame for db_feature_table if 'db = other'")
  }
  reads <- obtain_reads(input_file, input_type, aligner, blast_fastas, quiet)
  unmapped <- is.na(reads[[1]]$rname)
  if (db == "ncbi") reads[[1]]$rname <- identify_rnames(reads)
  mapped_rname <- as.character(reads[[1]]$rname[!unmapped])
  mapped_qname <- reads[[1]]$qname[!unmapped]
  mapped_cigar <- reads[[1]]$cigar[!unmapped]
  mapped_qwidth <- reads[[1]]$qwidth[!unmapped]
  if (blast_fastas) mapped_seqs <- reads[[1]]$seq[!unmapped]
  if (aligner == "bowtie2") {
    # mapped alignments used
    map_edit_or_align <- reads[[1]][["tag"]][["AS"]][!unmapped]
  } else if (aligner == "subread") {
    map_edit_or_align <-
      reads[[1]][["tag"]][["NM"]][!unmapped] # mapped edits used
  }
  read_names <- unique(mapped_qname)
  accessions <- as.character(unique(mapped_rname))
  if (db == "ncbi") {
    tax_id_all <- find_accessions(accessions, NCBI_key, quiet = quiet)
    taxids <- vapply(tax_id_all, function(x) x[1], character(1)) %>%
      missing_acc_finder(NCBI_key)
    genome_names <- vapply(tax_id_all, function(x) attr(x, "name"),
                           character(1))
    # Accession ids for any unknown genomes (likely removed from db)
    unk_inds <- which(is.na(taxids))
    genome_names[unk_inds] <- paste("unknown genome; accession ID is",
                                    accessions[unk_inds])
    taxids[unk_inds] <- accessions[unk_inds]
  } else if (db == "silva") {
    tax_id_all <- stringr::str_split(accessions, ";", n =2)
    taxids <- sapply(tax_id_all, `[[`, 1)
    genome_names <- sapply(tax_id_all, `[[`, 2)
    # Fix names
    mapped_rname <- stringr::str_split(mapped_rname, ";", n = 2) %>%
      sapply(`[[`, 1)
    accessions <- as.character(unique(mapped_rname))
  } else if (db == "other") {
    tax_id_all <- dplyr::tibble(`Feature ID` = accessions) %>%
      dplyr::left_join(db_feature_table, by = "Feature ID")
    taxids <- tax_id_all %>% dplyr::select(1) %>% unlist() %>% unname()
    genome_names <- tax_id_all %>% dplyr::select(2) %>% unlist() %>%
      unname()
  }
  unique_taxids <- unique(taxids)
  taxid_inds <- match(taxids, unique_taxids)
  unique_genome_names <- genome_names[!duplicated(taxid_inds)]
  if (!quiet) message("\tFound ", length(unique_taxids),
                      " unique taxa")
  # Make an aligment matrix (rows: reads, cols: unique taxids)
  if (!quiet) message("Setting up the EM algorithm")
  qname_inds <- match(mapped_qname, read_names)
  rname_inds <- match(mapped_rname, accessions)
  rname_tax_inds <- taxid_inds[rname_inds] #accession to taxid
  # Order based on read names
  rname_tax_inds <- rname_tax_inds[order(qname_inds)]
  cigar_strings <- mapped_cigar[order(qname_inds)]
  qwidths <- mapped_qwidth[order(qname_inds)]
  if (aligner == "bowtie2") {
    # mapped alignments used
    scores <- map_edit_or_align[order(qname_inds)]
  } else if (aligner == "subread") {
    # mapped edits used
    scores <- map_edit_or_align[order(qname_inds)]
  } else if (aligner == "other") scores <- 1
  qname_inds <- sort(qname_inds)
  exp_alignment_scores <- get_alignscore(aligner, cigar_strings,
                                         count_matches, scores, qwidths)
  combined <- dplyr::bind_cols("qname" = qname_inds,
                               "rname" = rname_tax_inds,
                               "scores" = exp_alignment_scores)
  results <- get_assignments(combined, convEM, maxitsEM, unique_taxids,
                             unique_genome_names, quiet = quiet)
  metascope_id_file <- results[[1]] %>% select("TaxonomyID", "Genome",
                                               "read_count", "Proportion",
                                               "readsEM", "EMProportion")
  utils::write.csv(metascope_id_file, file = out_file, row.names = FALSE)
  if (!quiet) message("Results written to ", out_file)

  if (blast_fastas){
    combined_distinct <- results[[2]]
    num_genomes <- min(num_genomes, nrow(results[[1]]))
    dir.create(file.path(out_dir, "fastas"))
    for (i in seq.int(1, num_genomes)) {
      current_genome <- results[[1]]$Genome[i] %>% stringr::word(1:2) %>% paste0(collapse = "_")
      current_rname_ind <- results[[1]]$hits_ind[i]
      read_indices <- combined_distinct %>% filter(rname == current_rname_ind) %>%
        pull(index)
      current_num_reads <- min(num_reads, length(read_indices))
      read_indices <- read_indices %>% sample(current_num_reads)
      seqs <- reads[[1]]$seq[read_indices]
      Biostrings::writeXStringSet(seqs, file.path(out_dir, "fastas",
                                                  paste0(i, "_", current_genome, ".fa")))
    }
  }

  if (update_bam){
    combined_single <- results[3]
    filter_which <- combined_single$single_hit
    bam_out <- file.path(out_dir, paste0(out_base, ".updated.bam"))
    Rsamtools::filterBam(file = input_file, destination = bam_out, filter = filter_which)
  }
  # Plotting genome locations
  num_plot <- num_species_plot
  if (is.null(num_plot)) num_plot <- min(nrow(results), 10)
  if (num_plot > 0) {
    if (!quiet) message("Creating coverage plots at ",
                        out_base, "_cov_plots")
    lapply(seq_along(results$TaxonomyID)[seq_len(num_plot)], function(x) {
      locations(as.numeric(results$TaxonomyID)[x],
                which_genome = results$Genome[x],
                accessions, taxids, reads, out_base, out_dir)})
  } else if (!quiet) message("No coverage plots created")
  return(results)
}

