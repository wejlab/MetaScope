globalVariables(c("qname", "rname"))

obtain_reads <- function(input_file, input_type, aligner){
    if (input_type == "bam") {
        message("Reading .bam file: ", input_file)
        to_pull <- c("qname", "rname", "cigar","qwidth", "pos")
        if (identical(aligner, "bowtie")) {
            params <- Rsamtools::ScanBamParam(what = to_pull, tag = c("AS"))
        } else if (identical(aligner, "subread")) {
            params <- Rsamtools::ScanBamParam(what = to_pull, tag = c("NM"))
        } else if (identical(aligner, "other")) {
            params <- Rsamtools::ScanBamParam(what = to_pull)
        }
        reads <- Rsamtools::scanBam(input_file, param = params)
    } else if (input_type == "rds") {
        reads <- list(as.list(readRDS(input_file)))
        if (identical(aligner, "bowtie")) {
            reads[[1]]$tag <- list("AS" = reads[[1]]$AS)
            reads[[1]]$AS <- NULL
        } else if (identical(aligner, "subread")) {
            reads[[1]]$tag <- list("NM" = reads[[1]]$NM)
            reads[[1]]$NM <- NULL
        }
    }
    return(reads)
}

find_accessions <- function(accessions, NCBI_key) {
    # Convert accessions to taxids and get genome names
    message("Obtaining taxonomy and genome names")
    # If URI length is greater than 2500 characters then split accession list
    URI_length <- nchar(paste(accessions, collapse = "+"))
    if (URI_length > 2500){
        chunks <- split(accessions, ceiling(seq_along(accessions) / 100))
        tax_id_all <- c()
        message(paste("Accession list broken into", length(chunks), "chunks"))
        for (i in 1:length(chunks)) {
            success <- FALSE
            attempt <- 0
            # Attempt to get taxid up to three times for each chunk
            while (!success) {
                try({
                    attempt <- attempt + 1
                    if (attempt > 1) message("Attempt #", attempt, " Chunk #", i)
                    suppressMessages(
                        tax_id_chunk <- taxize::genbank2uid(id = chunks[[i]],
                                                            key = NCBI_key))
                    Sys.sleep(1)
                    tax_id_all <- c(tax_id_all, tax_id_chunk)
                    success <- TRUE
                })
            }
        }
    } else suppressMessages(tax_id_all <- taxize::genbank2uid(id = accessions,
                                                              key = NCBI_key))
    return(tax_id_all)
}


#' Count the number of base lengths in a CIGAR string for a given operation
#' 
#' The ‘CIGAR’ (Compact Idiosyncratic Gapped Alignment Report) string is how
#' the SAM/BAM format represents spliced alignments. This function will accept
#' a CIGAR string for a single read and a single character indicating the
#' operation to be parsed in the string. An operation is a type of column that
#' appears in the alignment, e.g. a match or gap. The integer following the
#' operator specifies a number of consecutive operations. The `count_matches()`
#' will identify all occurrences of the operator in the string input, add them,
#' and return an integer number representing the total number of operations
#' for the read that was summarized by the input CIGAR string.
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
#' sapply(cigar3, count_matches, char = "I")
#'

count_matches <- function(x, char = "M") {
    if (length(char) != 1) {
        stop("Please provide a single character ",
             "operator with which to parse.")
    } else if (length(x) != 1) {
        stop("Please provide a single CIGAR string to be parsed.")
    }
    pattern <- paste("\\d+", char , sep = "")
    ind <- gregexpr(pattern, x)[[1]]
    start <- as.numeric(ind)
    end <- start + attr(ind, "match.length") - 2
    out <- sum(as.numeric(apply(cbind(start, end), 1,
                                function(y) substr(x, start = y[1],
                                                   stop = y[2]))))
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
#' @param input_file The path to the input file
#' 
#' @return A plot of the read coverage for a given genome

locations <- function(which_taxid, which_genome,
                      accessions, taxids, reads, input_file) {
    plots_save <- paste(stringr::str_split(tools::file_path_sans_ext(input_file),
                                           "\\.")[[1]][1],
                        "coverage", sep = "_")
    # map back to accessions
    choose_acc <- paste(accessions[which(as.numeric(taxids) %in% which_taxid)])
    # map back to BAM
    map2bam_acc <- which(reads[[1]]$rname %in% choose_acc)
    # Split genome name to make digestible
    this_genome <- strsplit(which_genome, " ")[[1]][1:2]
    use_name <- paste(this_genome, collapse = " ")
    coverage <- round(mean(seq_len(338099) %in% unique(
        reads[[1]]$pos[map2bam_acc])), 3)
    # Plotting
    dfplot <- dplyr::tibble(x = reads[[1]]$pos[map2bam_acc])
    ggplot2::ggplot(dfplot, ggplot2::aes(x)) + 
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::labs(main = paste("Positions of reads mapped to", use_name),
                      xlab = "Leftmost position in genome",
                      ylab = "Read Count",
                      subtitle = paste("Coverage is", coverage),
                      caption = paste0("Accession Number: ", choose_acc))
    ggplot2::ggsave(paste0(plots_save, "/",
                           stringr::str_replace(use_name, " ", "_"), ".png"),
                    device = "png")
}

#' Identify which genomes are represented in a sample
#'
#' This function will read in a .bam or .rds file, annotate the taxonomy and genome
#' names, reduce the mapping ambiguity using a mixture model, and output a
#' .csv file with the results. Currently, it assumes that the genome
#' library/.bam files use NCBI accession names for reference names (rnames in
#' .bam file). 
#'
#' @param input_file The .bam or .rds file that needs to be identified.
#' @param input_type Extension of file input. Should be either "bam" or "rds".
#' Default is "bam".
#' @param aligner The aligner which was used to create the bam file. Default is 
#' "subread" but can also be set to "bowtie" or "other"
#' @param NCBI_key (character) NCBI Entrez API key. optional.
#' See taxize::use_entrez(). Due to the high number of
#' requests made to NCBI, this function will be less prone to errors
#' if you obtain an NCBI key.
#' You may enter the string as an input or set it as ENTREZ_KEY in .Renviron.
#' @param out_file The name of the .csv output file. Defaults to the input_file
#' basename plus ".metascope_id.csv".
#' @param EMconv The convergence parameter of the EM algorithm. Default set at
#' \code{1/10000}.
#' @param EMmaxIts The maximum number of EM iterations, regardless of whether
#' the EMconv is below the threshhold. Default set at \code{50}.
#' If set at \code{0}, the algorithm skips the EM step and summarizes the .bam
#' file 'as is'
#' @param num_species_plot The number of genome coverage plots to be saved.
#' Default is \code{NULL}, which saves coverage plots for the ten most
#' highly abundant species.
#' @param deprecated_ind Logical. If you aligned using an old bowtie index, the bam
#' qnames may be different, and the function will not run without setting
#' this parameter to \code{TRUE}.
#' 
#' @return
#' This function returns a .csv file with annotated read counts to genomes with
#' mapped reads. The function itself returns the output .csv file name.
#'
#' @export
#'
#' @examples
#' #### Align reads to reference library and then apply metascope_id()
#' 
#' ## Assuming filtered bam files already exist
#' 
#' ## Subread aligned bam file
#' 
#' ## Create object with path to filtered subread bam file 
#' bamPath <- system.file("extdata","subread_target.filtered.bam",
#' package = "MetaScope")
#' 
#' ## Run metascope id with the aligner option set to subread
#' metascope_id(input_file = bamPath, aligner = "subread")
#' 
#' ## Bowtie aligned bam file 
#' 
#' ## Create object with path to filtered subread bam file 
#' bamPath <- system.file("extdata","bowtie_target.filtered.bam",
#' package = "MetaScope")
#' 
#' ## Run metascope id with the aligner option set to bowtie
#' metascope_id(input_file = bamPath, aligner = "bowtie")
#' 
#' ## Different or unknown aligned bam file
#' 
#' ## Create object with path to unknown origin bam file 
#' bamPath <- system.file("extdata","subread_target.filtered.bam",
#' package = "MetaScope")
#' 
#' ## Run metascope id with the aligner option set to other 
#' metascope_id(input_file = bamPath, aligner = "other")
#' 

metascope_id <- function(input_file, input_type = "bam", aligner = "subread",
                         NCBI_key = NULL,
                         out_file = paste(tools::file_path_sans_ext(input_file),
                                          ".metascope_id.csv", sep = ""),
                         EMconv = 1/10000, EMmaxIts = 25,
                         num_species_plot = NULL,
                         deprecated_ind = FALSE) {

    # Check to make sure valid aligner is specified
    if (aligner != "bowtie" && aligner != "subread" && aligner != "other")
        stop("Please make sure aligner is set to either 'bowtie', 'subread',",
             " or 'other'")
    reads <- obtain_reads(input_file, input_type, aligner)
    unmapped <- is.na(reads[[1]]$rname)
    # Account for potential index issues
    if (deprecated_ind) {
        # 2015
        mapped_2015 <- reads[[1]]$rname[!unmapped] %>%
            stringr::str_split(., pattern = "ref\\|", n = 2) %>%
            sapply(., function(x) x[2]) %>%
            stringr::str_split(., pattern = "\\|", n = 2) %>%
            sapply(., function(x) x[1])
        # 2018
        mapped_2018 <- reads[[1]]$rname[!unmapped] %>%
            stringr::str_split(., pattern = "ion\\|", n = 2) %>%
            sapply(., function(x) x[2])
        # Select_best
        ind <- sum(is.na(mapped_2015)) < sum(is.na(mapped_2018))
        if (ind) {
            mapped_rname <- mapped_2015
        } else mapped_rname <- mapped_2018
    } else {
        mapped_rname <- reads[[1]]$rname[!unmapped]
    }
    mapped_qname <- reads[[1]]$qname[!unmapped]
    mapped_cigar <- reads[[1]]$cigar[!unmapped]
    mapped_qwidth <- reads[[1]]$qwidth[!unmapped]

    if (aligner == "bowtie") {
        mapped_alignment <- reads[[1]][["tag"]][["AS"]][!unmapped]
    } else if (aligner == "subread")
        mapped_edit <- reads[[1]][["tag"]][["NM"]][!unmapped]

    read_names <- unique(mapped_qname)
    accessions <- as.character(unique(mapped_rname))
    message("\tFound ", length(read_names), " reads aligned to ",
            length(accessions), " NCBI accessions")
    tax_id_all <- find_accessions(accessions, NCBI_key)
    taxids <- vapply(tax_id_all, function(x) x[1], character(1))
    unique_taxids <- unique(taxids)
    taxid_inds <- match(taxids, unique_taxids)
    genome_names <- vapply(tax_id_all, function(x) attr(x, "name"),
                           character(1))
    unique_genome_names <- genome_names[!duplicated(taxid_inds)]
    message("\tFound ", length(unique_taxids), " unique NCBI taxonomy IDs")

    # Make an aligment matrix (rows: reads, cols: unique taxids)
    message("Setting up the EM algorithm")
    qname_inds <- match(mapped_qname, read_names)
    rname_inds <- match(mapped_rname, accessions)
    rname_tax_inds <- taxid_inds[rname_inds] #accession to taxid

    # Order based on read names
    rname_tax_inds <- rname_tax_inds[order(qname_inds)]
    cigar_strings <- mapped_cigar[order(qname_inds)]
    qwidths <- mapped_qwidth[order(qname_inds)]

    if (aligner == "bowtie") {
        scores <- mapped_alignment[order(qname_inds)]
    } else if (aligner == "subread") {
        scores <- mapped_edit[order(qname_inds)]
    } else if (aligner == "other") scores <- 1
    qname_inds <- sort(qname_inds)

    #Subread alignment scores: CIGAR string matches - edit score
    if (identical(aligner, "subread")) {
        num_match <- unlist(vapply(cigar_strings, count_matches,
                                   USE.NAMES = FALSE, double(1)))
        alignment_scores <- num_match - scores
        scaling_factor <- 100.0 / max(alignment_scores)
        relative_alignment_scores <- alignment_scores - min(alignment_scores)
        exp_alignment_scores <- exp(relative_alignment_scores * scaling_factor)
    } else if (identical(aligner, "bowtie")) {
        # Bowtie2 alignment scores: AS value + read length (qwidths) 
        alignment_scores <- scores + qwidths
        scaling_factor <- 100.0 / max(alignment_scores)
        relative_alignment_scores <- alignment_scores - min(alignment_scores)
        exp_alignment_scores <- exp(relative_alignment_scores * scaling_factor)
    } else if (identical(aligner, "other")) {
        # Other alignment scores: No assumptions
        exp_alignment_scores <- 1
    }
    combined <- dplyr::bind_cols("qname" = qname_inds,
                                 "rname" = rname_tax_inds,
                                 "scores" = exp_alignment_scores)
    input_distinct <- dplyr::distinct(combined, qname, rname, .keep_all = TRUE)
    qname_inds_2 <- input_distinct$qname
    rname_tax_inds_2 <- input_distinct$rname
    scores_2 <- input_distinct$scores
    non_unique_read_ind <- unique(combined[[1]][(
        duplicated(input_distinct[,1]) | duplicated(input_distinct[, 1],
                                                    fromLast = TRUE))])
    # 1 if read is multimapping, 0 if read is unique
    y_ind_2 <- as.numeric(unique(input_distinct[[1]]) %in% non_unique_read_ind)
    gammas <- Matrix::sparseMatrix(qname_inds_2, rname_tax_inds_2, x = scores_2)
    pi_old <- rep(1 / nrow(gammas), ncol(gammas))
    pi_new <- theta_new <- Matrix::colMeans(gammas)
    conv <- max(abs(pi_new - pi_old) / pi_old)
    it <- 0

    message("Starting EM iterations")
    while (conv > EMconv & it < EMmaxIts) {
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
        print(c(it, conv))
    }
    message("\tDONE! Converged in ", it, " iterations.")

    # Collect results
    hit_which <- qlcMatrix::rowMax(gammas_new, which = TRUE)$which
    best_hit <- Matrix::colSums(hit_which)
    names(best_hit) <- seq_along(best_hit)
    best_hit <- best_hit[best_hit != 0]
    hits_ind <- as.numeric(names(best_hit))
    final_taxids <- unique_taxids[hits_ind]
    final_genomes <- unique_genome_names[hits_ind]
    proportion <- best_hit / sum(best_hit)
    gammasums <- Matrix::colSums(gammas_new)
    EMreads <- round(gammasums[hits_ind], 1)
    EMprop <- gammasums[hits_ind] / sum(gammas_new)
    results <- cbind(TaxonomyID = final_taxids, Genome = final_genomes,
                     read_count = best_hit, Proportion = proportion,
                     EMreads = EMreads, EMProportion = EMprop)
    results <- as.data.frame(results[order(best_hit, decreasing = TRUE), ])
    message("Found reads for ", length(best_hit), " genomes")

    # PLotting of genome locations
    if (is.null(num_species_plot)){
        num_species_plot <- min(length(final_taxids), 10)
    }
    if (num_species_plot > 0) {
        message("Creating coverage plots")
        sapply(seq_along(results$TaxonomyID)[1:num_species_plot],
               function(x) locations(as.numeric(results$TaxonomyID)[x],
                                     which_genome = results$Genome[x],
                                     accessions, taxids, reads,
                                     input_file))
    }
    utils::write.csv(results, file = out_file, row.names = FALSE)
    message("Results written to ", out_file)
    return(results)
}
