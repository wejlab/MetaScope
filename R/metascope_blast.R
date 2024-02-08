library(Rsamtools)
library(taxize)
library(Biostrings)
library(dplyr)
library(R.utils)
library(rBLAST)
library(stringr)

#' Gets sequences from bam file
#'
#' Returns fasta sequences from a bam file with a given taxonomy ID
#'
#' @param id Taxonomy ID of genome to get sequences from
#' @param bam_file A sorted bam file and index file, loaded with Rsamtools::bamFile
#' @param n Number of sequences to retrieve
#'
#' @return Biostrings format sequences

get_seqs <- function(id, bam_file, n = 10) {
  # Get sequence info (Genome Name) from bam file
  seq_info_df <- data.frame(Rsamtools::seqinfo(bam_file))
  seq_info_df$seqnames <- row.names(seq_info_df)
  allGenomes <- grep(paste0("ti|", id), seq_info_df$seqnames, value = TRUE, fixed = TRUE)
  # Sample one of the Genomes that match
  Genome = sample(allGenomes, 1)
  # Scan Bam file for all sequences that match genome
  param <- Rsamtools::ScanBamParam(what = c("rname", "seq"),
                                   which = GRanges(Genome, IRanges(1, 1e+07)))
  allseqs <- Rsamtools::scanBam(bam_file, param = param)[[1]]
  n = min(n, length(allseqs$seq))
  seqs <- sample(allseqs$seq, n)

  return(seqs)
}


#' Converts NCBI taxonomy ID to scientific name
#'
#' @param taxids List of NCBI taxids to convert to scientific name
#' @param batch_size Number of taxids to submit to enterez (Need Key), max 10
#'
#' NCBI_key <- ""
#' options("ENTREZ_KEY" = NCBI_key)
#'
#' @return Returns a dataframe of blast results for a metascope result

taxid_to_name <- function(taxids, batch_size = 10) {
  nums <- length(taxids)
  i = 1
  df = data.frame(staxids = NULL, name = NULL)
  while (i <= length(taxids)) {
    tmp_df <- taxize::id2name(taxids[i:min(i+batch_size-1, nums)], db = "ncbi") %>%
      do.call(rbind, .) %>% rename(staxids = id) %>% select(staxids, name)
    df <- rbind(df, tmp_df)
    i = i + batch_size
    Sys.sleep(1)
  }
  df <- transform(df, staxids = as.integer(staxids))
  return(df)
}

#' rBLAST_single_result
#'
#' @param results_table A dataframe of the Metascope results
#' @param bam_file A sorted bam file and index file, loaded with Rsamtools::bamFile
#' @param which_result Index in results_table for which result to Blast search
#' @param num_reads Number of reads to blast per result
#' @param hit_list Number of how many blast results to fetch per read
#' @param num_threads Number of threads if multithreading
#' @param db_path Blast database path
#'
#' @return Returns a dataframe of blast results for a metascope result

rBLAST_single_result <- function(results_table, bam_file, which_result = 1, num_reads = 100,
                                 hit_list = 10, num_threads = 1, db_path, quiet = FALSE) {
  res <- tryCatch( #If any errors, should just skip the organism
    {
      genome_name <- results_table[which_result,2]
      if (!quiet) message("Current id: ", genome_name)
      tax_id <- results_table[which_result,1]
      if (!quiet) message("Current ti: ", tax_id)
      fasta_seqs <- get_seqs(id = tax_id, bam_file = bam_file, n = num_reads)
      blast_db <- rBLAST::blast(db = db_path, type = "blastn")
      blast_res <- rBLAST::predict(blast_db, fasta_seqs,
                                   custom_format ="qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids",
                                   BLAST_args = paste0("-max_target_seqs ", hit_list, " -num_threads ", num_threads))
      taxize_genome_df <- taxid_to_name(unique(blast_res$staxids))
      blast_res$MetaScope_Taxid <- tax_id
      blast_res$MetaScope_Genome <- genome_name
      blast_res <- dplyr::left_join(blast_res, taxize_genome_df, by = "staxids")
      blast_res
    },
    error = function(e) {
      cat("Error", conditionMessage(e), "\n")
      blast_res <- data.frame(qseqid=NA, sseqid=NA, pident=NA, length=NA,
                              mismatch=NA, gapopen=NA, qstart=NA, qend=NA,
                              sstart=NA, send=NA, evalue=NA, bitscore=NA, staxids=NA)
      tax_id <- results_table[which_result,1]
      genome_name <- results_table[which_result,2]
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
#' @param results_table A dataframe of the Metascope results
#' @param bam_file A sorted bam file and index file, loaded with Rsamtools::bamFile
#' @param num_results A number indicating number of Metascope results to blast
#' @param num_reads_per_result A number indicating number of reads to blast per result
#' @param hit_list A number of how many blast results to fetch for each read
#' @param num_threads Number of threads if multithreading
#' @param db_path Blast database path
#' @param out_path Output directory to save csv files, including base name of files
#'
#' @return Creates and exports num_results number of csv files with blast results from local blast

rBlast_results <- function(results_table, bam_file, num_results = 10, num_reads_per_result = 100, hit_list = 10,
                           num_threads = 1, db_path, out_path, sample_name = NULL) {
  for (i in seq.int(num_results)) {
    df <- rBLAST_single_result(results_table, bam_file, which_result = i,
                               num_reads = num_reads_per_result, hit_list = hit_list,
                               num_threads = num_threads, db_path = db_path)
    tax_id <- results_table[i,1]
    write.csv(df, file.path(out_path, paste0(sprintf("%05d", i), "_", sample_name, "_", "tax_id_", tax_id, ".csv")))
  }
}



#' Calculates result metrics from a blast results table
#'
#' This function calculates the best hit (genome with most blast read hits),
#' uniqueness score (total number of genomes hit), species percentage hit
#' (percentage of reads where MetaScope species also matched the blast hit
#' species), genus percentage hit (percentage of reads where blast genus matched
#' MetaScope aligned genus) and species contaminant score (percentage of reads
#' that blasted to other species genomes) and genus contaminant score (percentage
#' of reads that blasted to other genus genomes)
#'
#' @param blast_results_table_path path for blast results csv file
#'
#' @return a dataframe with best_hit, uniqueness_score, species_percentage_hit
#' genus_percentage_hit, species_contaminant_score, and genus_contaminant_score
#'

blast_result_metrics <- function(blast_results_table_path){
  tryCatch(
    {
      # Load in blast results table
      blast_results_table <- read.csv(blast_results_table_path, header = TRUE)

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
      blast_results_table <- blast_results_table %>%
        dplyr::mutate(.data$MetaScope_genus = word(.data$MetaScope_Genome, 1, 1, sep = " ")) %>%
        dplyr::mutate(.data$MetaScope_species = word(.data$MetaScope_Genome, 1, 2, sep = " ")) %>%
        dplyr::mutate(.data$query_genus = word(.data$name, 1, 1, sep = " ")) %>%
        dplyr::mutate(.data$query_species = word(.data$name, 1, 2, sep = " "))

      # Removing duplicate query num and query species
      blast_results_table <- blast_results_table %>%
        dplyr::distinct(.data$qseqid, .data$query_species, .keep_all = TRUE)

      # Calculating Metrics
      best_hit <- blast_results_table %>%
        dplyr::group_by(.data$query_species) %>%
        dplyr::summarise(.data$num_reads = n()) %>%
        dplyr::slice_max(.data$num_reads, with_ties = FALSE)

      uniqueness_score <- blast_results_table %>%
        dplyr::group_by(.data$query_species) %>%
        dplyr::summarise(.data$num_reads = n()) %>%
        nrow()

      species_percentage_hit <- blast_results_table %>%
        dplyr::filter(.data$MetaScope_species == .data$query_species) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      genus_percentage_hit <- blast_results_table %>%
        dplyr::filter(.data$MetaScope_genus == .data$query_genus) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      species_contaminant_score <- blast_results_table %>%
        dplyr::filter(.data$MetaScope_species != .data$query_species) %>%
        dplyr::distinct(.data$qseqid, .keep_all = TRUE) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      genus_contaminant_score <- blast_results_table %>%
        dplyr::filter(.data$MetaScope_genus != .data$query_genus) %>%
        dplyr::distinct(.data$qseqid, .keep_all = TRUE) %>%
        nrow() / length(unique(blast_results_table$qseqid))

      return(data.frame(best_hit = best_hit$query_species,
                        uniqueness_score = uniqueness_score,
                        species_percentage_hit = species_percentage_hit,
                        genus_percentage_hit = genus_percentage_hit,
                        species_contaminant_score = species_contaminant_score,
                        genus_contaminant_score = genus_contaminant_score))
    },
    error = function(e)
    {
      cat("Error", conditionMessage(e), "/n")
      return(data.frame(best_hit = 0,
                        uniqueness_score = 0,
                        species_percentage_hit = 0,
                        genus_percentage_hit = 0,
                        species_contaminant_score = 0,
                        genus_contaminant_score = 0))
    }
  )
}




#' Blast reads from MetaScope aligned files
#'
#' This is the main MetaScope Blast function. It aligns the top MetaScope
#' results to NCBI BLAST database. It requires that BLAST and a separate
#' nucleotide database is installed. This function BLAST reads from the MetaScope
#' aligned reads and returns a csv file updated with BLAST result metrics.
#'
#' @param metascope_id_path Path of MetaScope id csv file
#' @param tmp_dir Path for Bam directory
#' @param out_dir Path for output directory
#' @param sample_name Sample name for output files
#' @param num_reads Max number of reads to blast per microbe
#' @param hit_list Character for number of blast hit results to keep
#' @param num_threads Number of threads if multithreading
#' @param db_path Blast database path
#'
#' @returns This function writes an updated csv file with metrics evaluating
#'
#' @export
#'
#' @examples
#' \donttest{
#' ## Run metascope id
#' ### Create temporary directory
#' file_temp <- tempfile()
#' dir.create(file_temp)
#'
#' bamPath <- system.file("extdata", "bowtie_target.bam", package = "MetaScope")
#' file.copy(bamPath, file_temp)
#'
#' metascope_id(file.path(file_temp, bamPath), aligner = "bowtie2",
#'   input_type = "bam")
#'
#' ## Run metascope blast
#' ### Get export name and metascope id results
#' out_base <- bamPath %>% base::basename() %>% strsplit(split = "\\.") %>%
#'     magrittr::extract2(1) %>% magrittr::extract(1)
#' metascope_id_path <- file.path(out_dir, paste0(out_base, ".metascope_id.csv"))
#'
#' metascope_blast(metascope_id_path, tmp_dir = file_tmp,
#'   out_dir = dirname(input_file), sample_name = out_base, db_path = ????)
#' }
#'
#' ## Remove temporary directory
#' unlink(file_temp, recursive = TRUE)
#'

metascope_blast <- function(metascope_id_path, ref_dir, out_dir, sample_name,
                            num_results = 10, num_reads = 100, hit_list = 10,
                            num_threads = 1, db_path) {

  # Sort and index bam file
  bam_file_path <- list.files(path = tmp_dir, pattern = "\\.bam$", full.names = TRUE)
  sorted_bam_file_path <- file.path(tmp_dir, paste0(sample_name, "_sorted"))
  Rsamtools::sortBam(bam_file_path, destination = sorted_bam_file_path)
  sorted_bam_file <- paste0(sorted_bam_file_path, ".bam")
  Rsamtools::indexBam(sorted_bam_file)
  bam_file <- Rsamtools::BamFile(sorted_bam_file, index = sorted_bam_file)

  # Load in metascope id file and clean unknown genomes
  metascope_id <- read.csv(metascope_id_path, header = TRUE)

  # Create blast directory in tmp directory to save blast results in
  blast_tmp_dir <- file.path(tmp_dir,"blast")
  dir.create(blast_tmp_dir)

  # Run rBlast on all metascope microbes
  rBlast_results(results_table = metascope_id, bam_file = bam_file, num_results = num_results,
                 num_reads_per_result = num_reads, hit_list = hit_list, num_threads = num_threads,
                 db_path = db_path, out_path = blast_tmp_dir, sample_name = sample_name)

  # Run Blast metrics
  blast_result_metrics_list <- lapply(list.files(blast_tmp_dir, full.names = TRUE),
                                      blast_result_metrics)

  # Append Blast Metrics to MetaScope results
  blast_result_metrics_df <- as.data.frame(do.call(rbind, blast_result_metrics_list))
  blast_result_metrics_df[(nrow(blast_result_metrics_df)+1):nrow(metascope_id),] <- NA

  metascope_blast_df <- data.frame(metascope_id, blast_result_metrics_df)
  write.csv(metascope_blast_df, file.path(out_dir, paste0(sample_name, ".metascope_blast.csv")))

}
