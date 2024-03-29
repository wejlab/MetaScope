#' Gets sequences from bam file
#'
#' Returns fasta sequences from a bam file with a given taxonomy ID
#'
#' @param id Taxonomy ID of genome to get sequences from
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param n Number of sequences to retrieve
#' @param bam_seqs A list of the sequence IDs from the bam file
#' @inheritParams metascope_blast
#'
#' @return Biostrings format sequences

get_seqs <- function(id, bam_file, n = 10, quiet, NCBI_key = NULL,
                     bam_seqs) {
  rlang::is_installed("GenomicRanges")
  rlang::is_installed("IRanges")
  allGenomes <- stringr::str_detect(id, bam_seqs)
  # Get sequence info (Genome Name) from bam file
  seq_info_df <- as.data.frame(Rsamtools::seqinfo(bam_file)) |>
    tibble::rownames_to_column("seqnames")
  # Sample one of the Genomes that match
  Genome <- sample(seq_info_df$seqnames[allGenomes], 1)
  # Scan Bam file for all sequences that match genome
  param <- Rsamtools::ScanBamParam(what = c("rname", "seq"),
                                   which = GenomicRanges::GRanges(
                                     seq_info_df$seqnames[allGenomes],
                                     IRanges::IRanges(1, 1e+07)))
  allseqs <- Rsamtools::scanBam(bam_file, param = param)[[1]]
  n <- min(n, length(allseqs$seq))
  seqs <- sample(allseqs$seq, n)
  return(seqs)
}

#' Converts NCBI taxonomy ID to scientific name
#'
#' @param taxids List of NCBI taxids to convert to scientific name
#' @inheritParams metascope_blast
#' @importFrom rlang .data
#' @return Returns a dataframe of blast results for a metascope result

taxid_to_name <- function(taxids, NCBI_key = NULL) {
  smol_func <- function(x) {
    out <- class_taxon(x, NCBI_key, 5)
    if ("genus" %in% out$rank) {
      genus <- dplyr::filter(out, rank == "genus") |> dplyr::pull("name")
    } else genus <- NA
    if ("species" %in% out$rank) {
      species <- dplyr::filter(out, rank == "species") |> dplyr::pull("name")
    } else species <- NA
    tibble::tibble("staxids" = as.integer(x), genus, species)
  }
  all_combined <- plyr::adply(taxids, 1, smol_func, .id = NULL)
  return(all_combined)
}

#' rBLAST_single_result
#'
#' @param results_table A dataframe of the Metascope results
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param which_result Index in results_table for which result to Blast search
#' @param num_reads Number of reads to blast per result
#' @param hit_list Number of how many blast results to fetch per read
#' @param num_threads Number of threads if multithreading
#' @param db_path Blast database path
#' @param bam_seqs A list of the sequence IDs from the bam file
#' @inheritParams metascope_blast
#'
#' @return Returns a dataframe of blast results for a metascope result

rBLAST_single_result <- function(results_table, bam_file, which_result,
                                 num_reads = 100, hit_list = 10, 
                                 num_threads = 1, db_path, quiet = quiet,
                                 NCBI_key = NULL, bam_seqs) {
  res <- tryCatch({ #If any errors, should just skip the organism
    rlang::is_installed("rBLAST")
    genome_name <- results_table[which_result, 2]
    if (!quiet) message("Current id: ", genome_name)
    tax_id <- results_table[which_result, 1]
    if (!quiet) message("Current ti: ", tax_id)
    fasta_seqs <- get_seqs(id = tax_id, bam_file = bam_file, n = num_reads,
                           quiet = quiet, NCBI_key = NCBI_key,
                           bam_seqs = bam_seqs)
    blast_db <- rBLAST::blast(db = db_path, type = "blastn")
    this_format <- paste("qseqid sseqid pident length mismatch gapopen",
                         "qstart qend sstart send evalue bitscore staxids")
    predict.BLAST <- utils::getFromNamespace("predict.BLAST", "rBLAST")
    blast_res <- predict.BLAST(blast_db, fasta_seqs,
                               custom_format = this_format,
                               BLAST_args = paste("-max_target_seqs",
                                                  hit_list, "-num_threads",
                                                  num_threads))
    taxize_genome_df <- taxid_to_name(unique(blast_res$staxids), NCBI_key = NCBI_key)
    blast_res$MetaScope_Taxid <- tax_id
    blast_res$MetaScope_Genome <- genome_name
    blast_res <- dplyr::left_join(blast_res, taxize_genome_df, by = "staxids")
    blast_res
  },
  error = function(e) {
    cat("Error", conditionMessage(e), "\n")
    all_colnames <- stringr::str_split(this_format, " ")[[1]]
    blast_res <- matrix(NA, ncol = length(all_colnames),
                        dimnames = list(c(), all_colnames)) |> as.data.frame()
    tax_id <- results_table[which_result, 1]
    genome_name <- results_table[which_result, 2]
    blast_res$MetaScope_Taxid <- tax_id
    blast_res$MetaScope_Genome <- genome_name
    blast_res$name <- NA
    blast_res
  }
  )
  return(res)
}

#' rBlast_results
#'
#' @param results_table A data.frame of the Metascope results
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param num_results A number indicating number of Metascope results to blast
#' @param num_reads_per_result A number indicating number of reads to blast per
#'   result
#' @param hit_list A number of how many blast results to fetch for each read
#' @param num_threads Number of threads if multithreading
#' @param db_path Blast database path
#' @param out_path Output directory to save csv files, including base name of
#'   files
#' @inheritParams metascope_blast
#'
#' @return Creates and exports num_results number of csv files with blast
#'   results from local blast

rBlast_results <- function(results_table, bam_file, num_results = 10,
                           num_reads_per_result = 100, hit_list = 10,
                           num_threads = 1, db_path, out_path,
                           sample_name = NULL, quiet = quiet,
                           NCBI_key = NULL) {
  # Grab all identifiers
  seq_info_df <- as.data.frame(Rsamtools::seqinfo(bam_file)) |>
    tibble::rownames_to_column("seqnames")
  bam_seqs <- find_accessions(seq_info_df$seqnames,
                              NCBI_key = NCBI_key, quiet = quiet) |>
    plyr::aaply(1, function(x) x[1])
  # Grab results
  num_results2 <- min(num_results, nrow(results_table))
  run_res <- function(i) {
    df <- rBLAST_single_result(results_table, bam_file, which_result = i,
                               num_reads = num_reads_per_result,
                               hit_list = hit_list, num_threads = num_threads,
                               db_path = db_path, quiet = quiet,
                               NCBI_key = NCBI_key, bam_seqs = bam_seqs)
    tax_id <- results_table[i, 1]
    utils::write.csv(df, file.path(out_path,
                                   paste0(sprintf("%05d", i), "_", sample_name,
                                          "_", "tax_id_", tax_id, ".csv")),
                     row.names = FALSE)
  }
  plyr::a_ply(seq_len(num_results2), 1, run_res)
}

#' Calculates result metrics from a blast results table
#'
#' This function calculates the best hit (genome with most blast read hits),
#' uniqueness score (total number of genomes hit), species percentage hit
#' (percentage of reads where MetaScope species also matched the blast hit
#' species), genus percentage hit (percentage of reads where blast genus matched
#' MetaScope aligned genus) and species contaminant score (percentage of reads
#' that blasted to other species genomes) and genus contaminant score
#' (percentage of reads that blasted to other genus genomes)
#'
#' @param blast_results_table_path path for blast results csv file
#' @inheritParams metascope_blast
#'
#' @return a dataframe with best_hit, uniqueness_score, species_percentage_hit
#'   genus_percentage_hit, species_contaminant_score, and
#'   genus_contaminant_score
#' 

blast_result_metrics <- function(blast_results_table_path, NCBI_key = NULL){
  tryCatch({
    # Load in blast results table
    blast_results_table <- utils::read.csv(blast_results_table_path, header = TRUE)
    
    # Remove any empty tables
    if (nrow(blast_results_table) < 2) {
      return(data.frame(best_hit = 0,
                        uniqueness_score = 0,
                        species_percentage_hit = 0,
                        genus_percentage_hit = 0,
                        species_contaminant_score = 0,
                        genus_contaminant_score = 0))
    }
    
    # Adding Species and Genus columns
    meta_tax <- taxid_to_name(unique(blast_results_table$MetaScope_Taxid),
                              NCBI_key = NCBI_key) |> dplyr::select(-"staxids")
    
    blast_results_table_2 <- blast_results_table |>
      dplyr::mutate("MetaScope_genus" = meta_tax$genus[1],
                    "MetaScope_species" = meta_tax$species[1]) |>
      dplyr::rename("query_genus" = "genus",
                    "query_species" = "species") |>
      # Getting best hit per read
      dplyr::group_by(.data$qseqid) |>
      dplyr::slice_max(.data$evalue, with_ties = TRUE) |>
      # Removing duplicate query num and query species
      dplyr::distinct(.data$qseqid, .data$query_species, .keep_all = TRUE)
    
    # Calculating Metrics
    best_hit <- blast_results_table_2 |>
      dplyr::group_by(.data$query_species) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      dplyr::slice_max(.data$num_reads, with_ties = FALSE)
    
    best_hit_strain <- blast_results_table_2 |>
      tidyr::separate("sseqid", c(NA, "gi", NA, NA, NA), sep = "\\|")  |>
      dplyr::group_by(.data$gi) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      dplyr::slice_max(.data$num_reads, with_ties = TRUE)
    if (nrow(best_hit_strain) > 1) {
      best_strain <- NA
    } else if (nrow(best_hit_strain) == 1) {
      res <- taxize::genbank2uid(best_hit_strain$gi, NCBI_key = NCBI_key)
      best_strain <- attr(res[[1]], "name")
    }
    
      all_results <- blast_results_table_2 |>
      dplyr::ungroup() |>
      tidyr::replace_na(replace = list("query_genus" = "Unknown",
                                       "query_species" = "Unknown",
                                       "MetaScope_genus" = "Unknown",
                                       "MetaScope_species" = "Unknown")) |>
      dplyr::mutate("is_equiv_sp" = .data$MetaScope_species == .data$query_species,
                    "is_equiv_gn" = .data$MetaScope_genus == .data$query_genus) |>
      dplyr::summarise(uniqueness_score = length(unique(.data$query_species)),
                       species_percentage_hit = mean(.data$is_equiv_sp),
                       genus_percentage_hit = mean(.data$is_equiv_gn)) |>
      dplyr::mutate(species_contaminant_score = 1 - .data$species_percentage_hit,
                    genus_contaminant_score = 1 - .data$genus_percentage_hit,
                    best_hit = best_hit$query_species,
                    best_hit_strain = best_strain)
    
    return(all_results)
  },
  error = function(e)
  {
    cat("Error", conditionMessage(e), "/n")
    return(data.frame(best_hit = 0,
                      uniqueness_score = 0,
                      species_percentage_hit = 0,
                      genus_percentage_hit = 0,
                      species_contaminant_score = 0,
                      genus_contaminant_score = 0,
                      best_hit = NA,
                      best_hit_strain = NA))
  }
  )
}

#' Blast reads from MetaScope aligned files
#'
#' This function allows the user to check a subset of identified reads against
#' NCBI BLAST and the nucleotide database to confirm or contradict results
#' provided by MetaScope. It aligns the top `metascope_id()` results to NCBI
#' BLAST database. It REQUIRES that command-line BLAST and a separate nucleotide
#' database have already been installed on the host machine. It returns a csv
#' file updated with BLAST result metrics.
#' 
#' This function assumes that you used the NCBI nucleotide database to process
#' samples, with a download date of 2021 or later. This is necessary for
#' compatibility with the bam file headers.
#'
#' This is highly computationally intensive and should be ran with multiple
#' cores, submitted as a multi-threaded computing job if possible.
#' 
#' Note, if best_hit_strain is FALSE, then no strain was observed more often
#' among the BLAST results.
#'
#' @param metascope_id_path Character string; path to a csv file output by
#'   `metascope_id()`.
#' @param bam_file_path Character string; full path to bam file for the same
#'   sample processed by `metascope_id`. Note that the `filter_bam` function
#'   must have output a bam file, and not a .csv.gz file. See
#'   `?filter_bam_bowtie` for more details. Defaults to
#'   \code{list.files(file_temp, ".updated.bam$")[1]}.
#' @param tmp_dir Character string, a temporary directory in which to host
#'   files.
#' @param out_dir Character string, path to output directory.
#' @param sample_name Character string, sample name for output files.
#' @param num_results Integer, the maximum number of taxa from the metascope_id
#'   output to check reads. Takes the top n taxa, i.e. those with largest
#'   abundance. Default is 10.
#' @param num_reads Integer, the maximum number of reads to blast per microbe.
#'   If the true number of reads assigned to a given taxon is smaller, then the
#'   smaller number will be chosen. Default is 100. Too many reads will involve
#'   more processing time.
#' @param hit_list Integer, number of blast hit results to keep. Default is 10.
#' @param num_threads Integer, number of threads if running in parallel
#'   (recommended). Default is 1.
#' @param db_path Character string. The database file to be searched (including
#'   basename, but without file extension). For example, if the nt database is
#'   in the nt folder, this would be /filepath/nt/nt assuming that the database
#'   files have the nt basename. Check this path if you get an error message
#'   stating "No alias or index file found".
#' @param quiet Logical, whether to print out more informative messages. Default
#'   is FALSE.
#' @param NCBI_key (character) NCBI Entrez API key. optional. See
#'   taxize::use_entrez(). Due to the high number of requests made to NCBI, this
#'   function will be less prone to errors if you obtain an NCBI key.
#'
#' @returns This function writes an updated csv file with metrics evaluating...
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ### Create temporary directory
#' file_temp <- tempfile()
#' dir.create(file_temp)
#' 
#' bamPath <- system.file("extdata", "bowtie_target.bam",
#'                        package = "MetaScope")
#' file.copy(bamPath, file_temp)
#' 
#' metascope_id(file.path(file_temp, "bowtie_target.bam"), aligner = "bowtie2",
#'              input_type = "bam", out_dir = file_temp, num_species_plot = 0,
#'              update_bam = TRUE)
#' 
#' ## Run metascope blast
#' ### Get export name and metascope id results
#' out_base <- bamPath %>% base::basename() %>% strsplit(split = "\\.") %>%
#'   magrittr::extract2(1) %>% magrittr::extract(1)
#' metascope_id_path <- file.path(file_temp, paste0(out_base, ".metascope_id.csv"))
#' 
#' # NOTE: change db_path to the location where your BLAST database is stored!
#' db <- "/restricted/projectnb/pathoscope/data/blastdb/nt/nt"
#' 
#' Sys.setenv(ENTREZ_KEY = "<your id here>")
#' 
#' metascope_blast(metascope_id_path,
#'                 bam_file_path = file.path(file_temp, "bowtie_target.updated.bam"),
#'                 tmp_dir = file_temp,
#'                 out_dir = file_temp,
#'                 sample_name = out_base,
#'                 db_path = db,
#'                 num_results = 10,
#'                 num_reads = 5,
#'                 hit_list = 10,
#'                 num_threads = 3,
#'                 quiet = FALSE,
#'                 NCBI_key = NULL)
#' 
#' ## Remove temporary directory
#' unlink(file_temp, recursive = TRUE)
#'}
#'

metascope_blast <- function(metascope_id_path,
                            bam_file_path = list.files(tmp_dir, ".updated.bam$",
                                                       full.names = TRUE)[1],
                            tmp_dir, out_dir, sample_name,
                            num_results = 10, num_reads = 100, hit_list = 10,
                            num_threads = 1, db_path, quiet = FALSE,
                            NCBI_key = NULL) {
  if (!is.numeric(threads)) threads <- 1
  # Sort and index bam file
  sorted_bam_file_path <- file.path(tmp_dir, paste0(sample_name, "_sorted"))
  Rsamtools::sortBam(bam_file_path, destination = sorted_bam_file_path)
  sorted_bam_file <- paste0(sorted_bam_file_path, ".bam")
  Rsamtools::indexBam(sorted_bam_file)
  bam_file <- Rsamtools::BamFile(sorted_bam_file, index = sorted_bam_file)
  
  # Load in metascope id file and clean unknown genomes
  metascope_id_in <- utils::read.csv(metascope_id_path, header = TRUE)
  
  # Create blast directory in tmp directory to save blast results in
  blast_tmp_dir <- file.path(tmp_dir, "blast")
  if(!dir.exists(blast_tmp_dir)) dir.create(blast_tmp_dir, recursive = TRUE)
  
  # Run rBlast on all metascope microbes
  rBlast_results(results_table = metascope_id_in, bam_file = bam_file,
                 num_results = num_results, num_reads_per_result = num_reads,
                 hit_list = hit_list, num_threads = num_threads,
                 db_path = db_path, out_path = blast_tmp_dir,
                 sample_name = sample_name, quiet = quiet, NCBI_key = NCBI_key)
  
  # Run Blast metrics
  blast_result_metrics_df <- plyr::adply(
    list.files(blast_tmp_dir, full.names = TRUE), 1, blast_result_metrics,
    NCBI_key = NCBI_key)
  
  # Append Blast Metrics to MetaScope results
  if (nrow(metascope_id_in) > nrow(blast_result_metrics_df)) {
    ind <- seq(nrow(blast_result_metrics_df) + 1, nrow(metascope_id_in))
    blast_result_metrics_df[ind, ] <- NA
  }
  print_file <- file.path(out_dir, paste0(sample_name, ".metascope_blast.csv"))
  metascope_blast_df <- data.frame(metascope_id_in, blast_result_metrics_df) |>
  # Add MetaScope original species after genome
    dplyr::mutate("MetaID Species" = taxid_to_name(metascope_id_in$TaxonomyID,
                                NCBI_key = NCBI_key)$species) |>
    dplyr::relocate("MetaID Species", .after = "Genome") |>
    dplyr::rename("MetaID Genome" = "Genome")
  utils::write.csv(metascope_blast_df, print_file)
  message("Results written to ", print_file)
  return(metascope_blast_df)
}