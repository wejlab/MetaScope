#' Adds in taxa if silva database

#' Returns MetaScope Table with silva taxa in separate columns

#' @param metascope_id_in MetaScope ID file with silva taxa
#' @param caching Boolean for if all_silva_headers.rds is already downloaded
#' @param path_to_write Path to save all_silva_headers.rds

add_in_taxa <- function(metascope_id_in, caching, path_to_write) {
  location <- "https://github.com/wejlab/metascope-docs/raw/main/all_silva_headers.rds"
  filename <- "all_silva_headers.rds"
  if (!caching) {
    if (!dir.exists(path_to_write)) dir.create(path_to_write)
    destination <- paste(path_to_write, filename, sep = "/")
    utils::download.file(location, destination)
  } else if (caching) {
    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, filename, "rname")$rid
    if (!length(rid)) {
      rid <- names(BiocFileCache::bfcadd(bfc, filename, location))
    }
    if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid))) {
      BiocFileCache::bfcdownload(bfc, rid)
      BiocFileCache::bfcrpath(bfc, rids = rid)
    } else {
      message("Caching is set to TRUE, ",
              "and it appears that this file is already downloaded ",
              "in the cache. It will not be downloaded again.")
    }
    destination <- BiocFileCache::bfcrpath(bfc, rids = rid)
  }
  all_silva_headers <- readRDS(destination)
  tax_table_pre <- metascope_id_in |>
    dplyr::mutate("TaxonomyID" = as.character(!!dplyr::sym("TaxonomyID"))) |>
    dplyr::distinct(!!dplyr::sym("TaxonomyID"), .keep_all = TRUE) |>
    dplyr::left_join(all_silva_headers, by = c("TaxonomyID")) |>
    dplyr::relocate("genus", "species", .after = "family") |>
    dplyr::relocate("read_count", "Proportion", "readsEM", "EMProportion", .after = "species") |>
    dplyr::select(-"Genome")
  return(tax_table_pre)
}

#' Adds in taxa if NCBI database

#' Returns MetaScope Table with NCBI taxa in separate columns

#' @param metascope_id_in MetaScope ID file with NCBI taxa qnames
#' @param NCBI_key (character) NCBI Entrez API key. optional.
#' See taxize::use_entrez(). Due to the high number of requests made to NCBI,
#' this function will be less prone to errors if you obtain an NCBI key.
#'

add_in_taxa_ncbi <- function(metascope_id_in, NCBI_key, param) {
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class",
                   "order", "family", "genus", "species")
  all_ncbi <- plyr::llply(metascope_id_in$TaxonomyID, .fun = class_taxon,
                          NCBI_key = NCBI_key, .progress = "text",
                          num_tries = 3)
  # fix unknowns
  unk_tab <- data.frame(name = "unknown", rank = taxon_ranks, id = 0)
  all_ncbi <- lapply(all_ncbi, function(x) if (all(is.na(x))) unk_tab else x)
  # Create table
  taxonomy_table <- BiocParallel::bplapply(all_ncbi, mk_table,
                                           taxon_ranks = taxon_ranks,
                                           BPPARAM = param) |>
    dplyr::bind_rows() |> as.data.frame() |>
    magrittr::set_colnames(taxon_ranks) |>
    dplyr::mutate(TaxonomyID = metascope_id_in$TaxonomyID) |>
    dplyr::relocate("TaxonomyID") |>
    dplyr::distinct(!!dplyr::sym("TaxonomyID"), .keep_all = TRUE) |>
    dplyr::left_join(metascope_id_in, by = c("TaxonomyID")) |>
    dplyr::relocate("genus", "species", .after = "family") |>
    dplyr::relocate("read_count", "Proportion", "readsEM", "EMProportion", .after = "species") |>
    dplyr::select(-"Genome")
  return(taxonomy_table)
}

#' Gets sequences from bam file
#'
#' Returns fasta sequences from a bam file with a given taxonomy ID
#'
#' @param id Taxonomy ID of genome to get sequences from
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param n Number of sequences to retrieve
#' @param bam_seqs A list of the sequence IDs from the bam file
#'
#' @return Biostrings format sequences

get_seqs <- function(id, bam_file, n = 10, bam_seqs) {
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


#' Gets multiple sequences from different accessions in a bam file
#'
#' Returns fasta sequences from a bam file with given taxonomy IDs
#'
#' @param ids_n List of vectors with Taxonomy IDs and the number of sequences to get from each
#' @param bam_file A sorted bam file and index file, loaded with
#'   Rsamtools::bamFile
#' @param seq_info_df Dataframe of sequence information from metascope_blast()
#' @param metascope_id_tax Data.frame of taxonomy information
#' @param sorted_bam_file Filepath to sorted bam file
#'
#' @return Biostrings format sequences
get_multi_seqs <- function(ids_n, bam_file, seq_info_df,
                           metascope_id_tax,
                           sorted_bam_file) {
  id <- ids_n[1]
  requireNamespace("GenomicRanges")
  requireNamespace("IRanges")
  allGenomes <- seq_info_df |>
    dplyr::filter(grepl(id, !!dplyr::sym("seqnames"), fixed = TRUE)) |>
    dplyr::pull("original_seqids")
  # Sample one of the Genomes that match
  Genome <- sample(allGenomes, 1)
  # Define BAM file with smaller yield size
  YS <- 1000
  bf <- Rsamtools::BamFile(sorted_bam_file, yieldSize = YS)
  
  alt_n <- metascope_id_tax |>
    dplyr::filter(!!dplyr::sym("TaxonomyID") == id) |> dplyr::pull("read_count")
  n_seqs <- ifelse(test = is.na(as.numeric(ids_n[[1]][2])),
                   yes = min(alt_n, 100),
                   no = as.numeric(ids_n[[1]][2]))
  
  YIELD <- function(bf) {
    param <- Rsamtools::ScanBamParam(what = c("rname", "seq"),
                                     which = GenomicRanges::GRanges(
                                       Genome,
                                       IRanges::IRanges(1, 1e+07)))
    value <- Rsamtools::scanBam(bf, param = param)[[1]]
    return(value)
  }
  MAP <- function(value) {
    if(length(value$rname) == 0) return(dplyr::tibble(NULL))
    y <- dplyr::tibble("UID" = value$rname |> unlist(),
                       "Sequence" = as.character(value$seq))
    return(y)
  }
  REDUCE <- function(x, y) {
    all_join <- dplyr::bind_rows(x, y)
    # Print all current reads to file
    return(all_join)
  }
  
  reduceByYield_multiseqs <- function(X, YIELD, MAP, REDUCE, n_end) {
    if (!Rsamtools::isOpen(X)) {
      Rsamtools::open.BamFile(X)
      on.exit(Rsamtools::close.BamFile(X))
    }
    result <- dplyr::tibble(NULL)
    
    while(nrow(result) < n_end) {
      data <- YIELD(X)
      result <- REDUCE(result, MAP(data))
      #print(nrow(result))
      #print(n_end)
    }
    return(result)
  }
  
  seqs <- reduceByYield_multiseqs(bf, YIELD, MAP, REDUCE, n_end = n_seqs)
  seqs <- seqs[sample(nrow(seqs), n_seqs), ]
  
  final_seqs <- list(rname = seqs$UID, seq = Biostrings::DNAStringSet(seqs$Sequence))
  return(final_seqs)
}

#' Converts NCBI taxonomy ID to scientific name
#'
#' @param taxids List of NCBI taxids to convert to scientific name
#' @param accessions_path Path to accessionsTaxa.sql
#' @importFrom rlang .data
#' @return Returns a dataframe of blast results for a metascope result
#'

taxid_to_name <- function(taxids, accessions_path) {
  taxids <- stringr::str_replace(taxids, ";(.*)$", "") |> as.integer()
  out <- taxonomizr::getTaxonomy(taxids, accessions_path)
  out_df <- as.data.frame(out) |> tibble::rownames_to_column("staxid") |>
    dplyr::select("staxid", "genus", "species") |>
    dplyr::mutate("species" = stringr::str_replace(!!dplyr::sym("species"), paste0(!!dplyr::sym("genus"), " "), "")) |>
    dplyr::mutate("staxid" = gsub("\\.", "", !!dplyr::sym("staxid"))) |>
    dplyr::mutate("staxid" = gsub("X", "", !!dplyr::sym("staxid"))) |>
    dplyr::mutate("staxid" = as.integer(!!dplyr::sym("staxid")))
  return(out_df)
}

#' Checks if blastn is installed
#'
#' @title Check if blastn exists on the system
#' @description This is an internal function that is not meant to be used
#'   outside of the package. It checks whether blastn exists on the system.
#' @return Returns TRUE if blastn exists on the system, else FALSE.
#'

check_blastn_exists <- function() {
  if (file.exists(Sys.which("blastn"))) return(TRUE)
  return(FALSE)
}

# BLASTn sequences


blastn_seqs <- function(db_path, fasta_path, res_path, hit_list, num_threads) {
  if (!check_blastn_exists()) {
    stop("blast not found, please install blast")
  }
  sys::exec_wait(
    "blastn", c("-db", db_path, "-query", fasta_path, "-out", res_path,
                "-outfmt", paste("10 qseqid sseqid pident length mismatch gapopen",
                                 "qstart qend sstart send evalue bitscore staxid"),
                "-max_target_seqs", hit_list, "-num_threads", num_threads,
                "-task", "megablast"))
}



#' blastn_single_result
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
#' @param out_path Path to output results.
#' @param fasta_dir Path to where fasta files are stored.
#' @inheritParams metascope_blast
#'
#' @return Returns a dataframe of blast results for a metascope result
#'

blastn_single_result <- function(results_table, bam_file, which_result,
                                 num_reads = 100, hit_list = 10,
                                 num_threads, db_path, quiet,
                                 accessions_path, bam_seqs, out_path,
                                 sample_name, fasta_dir = NULL) {
  res <- tryCatch({ #If any errors, should just skip the organism
    genome_name <- results_table[which_result, 8]
    if (!quiet) message("Current id: ", genome_name)
    tax_id <- results_table[which_result, 14] |> stringr::str_split(",") |>
      dplyr::first() |> dplyr::first() # Grabs First TaxID
    if (!quiet) message("Current ti: ", tax_id)
    # Generate sequences to blast
    if (!is.null(fasta_dir)) {
      fasta_path <- list.files(path = fasta_dir, full.names = TRUE)[which_result]
    } else {
      fasta_path <- get_seqs(id = tax_id, bam_file = bam_file, n = num_reads,
                             bam_seqs = bam_seqs)
    }
    
    res_path = file.path(out_path, paste0(sprintf("%05d", which_result), "_",
                                          sample_name, ".csv"))
    blastn_seqs(db_path, fasta_path, res_path = res_path, hit_list, num_threads)
    blast_res <- utils::read.csv(res_path, header = FALSE) |>
      magrittr::set_colnames(c("qseqid", "sseqid", "pident", "length", "mismatch",
                               "gapopen","qstart", "qend", "sstart", "send",
                               "evalue", "bitscore", "staxid")) |>
      dplyr::mutate("MetaScope_Taxid" = tax_id,
                    "MetaScope_Genome" = as.character(genome_name))
      taxize_genome_df <- taxid_to_name(unique(blast_res$staxid),
                                        accessions_path = accessions_path)
    blast_res <- dplyr::left_join(blast_res, taxize_genome_df, by = "staxid")
    blast_res
  },
  error = function(e) {
    cat("Error", conditionMessage(e), "\n")
    this_format <- paste("qseqid sseqid pident length mismatch gapopen",
                         "qstart qend sstart send evalue bitscore staxid")
    all_colnames <- stringr::str_split(this_format, " ")[[1]]
    blast_res <- matrix(NA, ncol = length(all_colnames),
                        dimnames = list(c(), all_colnames)) |> as.data.frame() |>
      dplyr::mutate("MetaScope_Taxid" = results_table[which_result, 1],
                    "MetaScope_Genome" =
                      as.character(results_table[which_result, 2]),
                    "name" = NA)
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
#' @param fasta_dir inc.
#' @inheritParams metascope_blast
#'
#' @return Creates and exports num_results number of csv files with blast
#'   results from local blast

blastn_results <- function(results_table, bam_file, num_results = 10,
                           num_reads_per_result = 100, hit_list = 10,
                           num_threads = 1, db_path, out_path,
                           sample_name = NULL, quiet = quiet,
                           accessions_path, fasta_dir = NULL,
                           NCBI_key = NULL, param) {
  # Grab all identifiers
  if (is.null(fasta_dir)) {
    seq_info_df <- as.data.frame(Rsamtools::seqinfo(bam_file)) |>
      tibble::rownames_to_column("seqnames")
    bam_seqs <- find_accessions(seq_info_df$seqnames,
                                NCBI_key = NCBI_key, quiet = quiet) |>
      plyr::aaply(1, function(x) x[1])
  } else {
    bam_seqs = NULL
  }
  # Grab results
  num_results2 <- min(num_results, nrow(results_table))
  run_res <- function(i) {
    df <- blastn_single_result(results_table, bam_file, which_result = i,
                               num_reads = num_reads_per_result,
                               hit_list = hit_list, num_threads = num_threads,
                               db_path = db_path, quiet = quiet, bam_seqs = bam_seqs,
                               out_path = out_path, sample_name = sample_name,
                               fasta_dir = fasta_dir, accessions_path = accessions_path)
    utils::write.csv(df, file.path(out_path,
                                   paste0(sprintf("%05d", i), "_", sample_name, ".csv")),
                     row.names = FALSE)
  }
  BiocParallel::bplapply(seq_len(num_results2), run_res, BPPARAM = param)
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
#' @return a vector with best_hit, uniqueness_score, species_percentage_hit
#'   genus_percentage_hit, species_contaminant_score, and
#'   genus_contaminant_score
#'

blast_result_metrics <- function(blast_results_table_path, accessions_path, db = NULL,
                                 NCBI_key = NULL){
  tryCatch({
    # Load in blast results table
    blast_results_table <- utils::read.csv(blast_results_table_path, header = TRUE)
    
    # Clean species results
    blast_results_table <- blast_results_table |>
      dplyr::filter(!grepl("sp.", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!grepl("uncultured", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!is.na(!!dplyr::sym("genus")))
    
    # Remove any empty tables
    if (nrow(blast_results_table) < 2) {
      error_vector <- c(rep(0, 5), rep(NA, 3)) |>
        magrittr::set_names(c("uniqueness_score", "species_percentage_hit", "genus_percentage_hit",
                            "species_contaminant_score", "genus_contaminant_score",
                            "best_hit_genus", "best_hit_species", "best_hit_strain"))
    return(error_vector)
    }
    # Clean species results
    blast_results_table <- blast_results_table |>
      dplyr::filter(!grepl("sp.", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!grepl("uncultured", !!dplyr::sym("species"), fixed = TRUE)) |>
      dplyr::filter(!is.na(!!dplyr::sym("genus")))
    # Adding MetaScope Species and Genus columns
    if (db == "silva") {
      # Cleaning up silva genome names
      #silva_genome <- gsub(";([a-z0-9])", " \\1", blast_results_table$MetaScope_Genome[1])
      silva_genome <- gsub(" uncultured", ";uncultured", blast_results_table$MetaScope_Genome[1])
      # Creating silva genus and species names
      silva_genus <- silva_genome |> strsplit(split = "_") |> sapply("[", 1)
      silva_species <- silva_genome |> strsplit(split = "_") |> sapply("[", 2)
      blast_results_table_2 <- blast_results_table |>
        dplyr::mutate("MetaScope_genus" = rep(silva_genus, nrow(blast_results_table)),
                      "MetaScope_species" = rep(silva_species, nrow(blast_results_table))) |>
        dplyr::rename("query_genus" = "genus",
                      "query_species" = "species") |>
        # Remove rows with NA
        tidyr::drop_na("query_genus", "query_species") |>
        # Getting best hit per read
        dplyr::group_by(!!dplyr::sym("qseqid")) |>
        dplyr::slice_min(!!dplyr::sym("evalue"), with_ties = TRUE) |>
        # Removing duplicate query num and query species
        dplyr::distinct(!!dplyr::sym("qseqid"), !!dplyr::sym("query_species"), .keep_all = TRUE)
      
    } else {
      meta_tax <- taxid_to_name(unique(blast_results_table$MetaScope_Taxid),
                                accessions_path = accessions_path) |>
        dplyr::select(-"staxid")
      blast_results_table_2 <- blast_results_table |>
        dplyr::mutate("MetaScope_genus" = meta_tax$genus[1],
                      "MetaScope_species" = meta_tax$species[1]) |>
        dplyr::rename("query_genus" = "genus",
                      "query_species" = "species") |>
        # Remove rows with NA
        tidyr::drop_na() |>
        # Getting best hit per read
        dplyr::group_by(!!dplyr::sym("qseqid")) |>
        dplyr::slice_min(!!dplyr::sym("evalue"), with_ties = TRUE) |>
        # Removing duplicate query num and query species
        dplyr::distinct(!!dplyr::sym("qseqid"), !!dplyr::sym("query_species"), .keep_all = TRUE)
    }
    # Calculating Metrics
    best_hit_genus <- blast_results_table_2 |>
      dplyr::group_by(!!dplyr::sym("query_genus")) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      dplyr::slice_max(!!dplyr::sym("num_reads"), with_ties = FALSE)
    
    best_hit_species <- blast_results_table_2 |>
      dplyr::group_by(!!dplyr::sym("query_species")) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      dplyr::slice_max(!!dplyr::sym("num_reads"), with_ties = FALSE)
    
    best_hit_strain <- blast_results_table_2 |>
      tidyr::separate("sseqid", c(NA, "gi", NA, NA, NA), sep = "\\|")  |>
      dplyr::group_by(!!dplyr::sym("gi")) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      dplyr::slice_max(!!dplyr::sym("num_reads"), with_ties = TRUE)
    if (nrow(best_hit_strain) > 1 | db == "silva") {
      best_strain <- NA
    } else if (nrow(best_hit_strain) == 1) {
      res <- taxize::genbank2uid(best_hit_strain$gi, NCBI_key = NCBI_key)
      best_strain <- attr(res[[1]], "name")
    }
    
    uniqueness_score <- blast_results_table_2 |>
      dplyr::group_by(!!dplyr::sym("query_species")) |>
      dplyr::summarise("num_reads" = dplyr::n()) |>
      nrow()
    
    species_percentage_hit <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_species") == !!dplyr::sym("query_species")) |>
      nrow() / length(unique(blast_results_table_2$qseqid))
    
    genus_percentage_hit <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_genus") == !!dplyr::sym("query_genus")) |>
      dplyr::distinct(!!dplyr::sym("qseqid"), .keep_all = FALSE) |>
      nrow() / length(unique(blast_results_table_2$qseqid))
    
    species_contaminant_score <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_species") != !!dplyr::sym("query_species")) |>
      dplyr::distinct(!!dplyr::sym("qseqid"), .keep_all = TRUE) |>
      nrow() / length(unique(blast_results_table_2$qseqid))
    
    genus_contaminant_score <- blast_results_table_2 |>
      dplyr::filter(!!dplyr::sym("MetaScope_genus") != !!dplyr::sym("query_genus")) |>
      dplyr::distinct(!!dplyr::sym("qseqid"), .keep_all = TRUE) |>
      nrow() / length(unique(blast_results_table_2$qseqid))
    
    
    all_results <- c(uniqueness_score, species_percentage_hit, genus_percentage_hit,
                     species_contaminant_score, genus_contaminant_score,
                     best_hit_genus$query_genus, best_hit_species$query_species,
                     best_strain)
    names(all_results) <- c("uniqueness_score", "species_percentage_hit", "genus_percentage_hit",
                            "species_contaminant_score", "genus_contaminant_score",
                            "best_hit_genus", "best_hit_species", "best_hit_strain")
    return(all_results)
  },
  error = function(e)
  {
    cat("Error", conditionMessage(e), "/n")
    error_vector <- c(0,0,0,0,0, NA, NA, NA)
    names(error_vector) <- c("uniqueness_score", "species_percentage_hit", "genus_percentage_hit",
                             "species_contaminant_score", "genus_contaminant_score",
                             "best_hit_genus", "best_hit_species", "best_hit_strain")
    
    return(error_vector)
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
#' @param db Currently accepts one of \code{c("ncbi", "silva", "other")}
#' Default is \code{"ncbi"}, appropriate for samples aligned against indices
#' compiled from NCBI whole genome databases. Alternatively, usage of an
#' alternate database (like Greengenes2) should be specified with
#' \code{"other"}.
#' @param fasta_dir Directory where fasta files for blast will be stored.
#' @param accessions_path Directory where accession files for blast are stored.
#'
#' @returns This function writes an updated csv file with metrics.
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
#' out_base <- bamPath |> base::basename() |> strsplit(split = "\\.") |>
#'   magrittr::extract2(1) |> magrittr::extract(1)
#' metascope_id_path <- file.path(file_temp, paste0(out_base, ".metascope_id.csv"))
#'
#' # NOTE: change db_path to the location where your BLAST database is stored!
#' db <- "/restricted/projectnb/pathoscope/data/blastdb/nt/nt"
#'
#' Sys.setenv(ENTREZ_KEY = "<your id here>")
#'
#'  metascope_blast(metascope_id_path,
#'               bam_file_path = file.path(file_temp, "bowtie_target.bam"),
#'               #bam_file_path = file.path(file_temp, "bowtie_target.updated.bam")
#'               tmp_dir = file_temp,
#'               out_dir = file_temp,
#'               sample_name = out_base,
#'               db_path = db_path,
#'               num_results = 10,
#'               num_reads = 5,
#'               hit_list = 10,
#'               num_threads = 3,
#'               db = "ncbi",
#'               quiet = FALSE,
#'               NCBI_key = NULL,
#'               fasta_dir = NULL ,
#'               accessions_path = NULL)
#'
#' ## Remove temporary directory
#' unlink(file_temp, recursive = TRUE)
#'}
#'

metascope_blast <- function(metascope_id_path,
                            bam_file_path = list.files(tmp_dir, ".updated.bam$",
                                                       full.names = TRUE)[1],
                            tmp_dir, out_dir, sample_name, fasta_dir = NULL,
                            num_results = 10, num_reads = 100, hit_list = 10,
                            num_threads = 1, db_path, quiet = FALSE,
                            NCBI_key = NULL, db = NULL, accessions_path = NULL) {
  # Parallelization settings
  if (!is.numeric(num_threads)) num_threads <- 1
  param <- BiocParallel::MulticoreParam(workers = num_threads)
  
  # Sort and index bam file
  sorted_bam_file_path <- file.path(tmp_dir, paste0(sample_name, "_sorted"))
  message("Sorting bam file to ", sorted_bam_file_path, ".bam")
  Rsamtools::sortBam(bam_file_path, destination = sorted_bam_file_path,
                     overwrite = TRUE)
  sorted_bam_file <- paste0(sorted_bam_file_path, ".bam")
  Rsamtools::indexBam(sorted_bam_file)
  bam_file <- Rsamtools::BamFile(sorted_bam_file, index = sorted_bam_file)
  
  # Load in metascope id file and clean unknown genomes
  metascope_id_in <- utils::read.csv(metascope_id_path, header = TRUE)
  
  # Group metascope id by species and create metascope species id
  if (db == "silva") {
    metascope_id_tax <- add_in_taxa(metascope_id_in, caching = FALSE,
                                    path_to_write = tmp_dir)
  } else if (db == "ncbi") {
    metascope_id_tax <- add_in_taxa_ncbi(metascope_id_in, param = param)
  }
  
  metascope_id_species <- metascope_id_tax |> dplyr::mutate("id" = dplyr::row_number()) |>
    dplyr::group_by(!!dplyr::sym("superkingdom"), !!dplyr::sym("kingdom"), !!dplyr::sym("phylum"),
                    !!dplyr::sym("class"),
                    !!dplyr::sym("order"), !!dplyr::sym("family"), !!dplyr::sym("genus"), !!dplyr::sym("species")) |>
    dplyr::summarise("read_count" = sum(!!dplyr::sym("read_count")),
                     "Proportion" = sum(!!dplyr::sym("Proportion")),
                     "readsEM" = sum(!!dplyr::sym("readsEM")),
                     "EMProportion" = sum(!!dplyr::sym("EMProportion")),
                     "IDs" = paste0(!!dplyr::sym("id"), collapse = ","),
                     "TaxonomyIDs" = paste0(!!dplyr::sym("TaxonomyID"), collapse = ","),
                     "read_proportions" = paste0(!!dplyr::sym("read_count")/sum(!!dplyr::sym("read_count")),
                                                 collapse = ","),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(!!dplyr::sym("read_count")))
  
  utils::write.csv(metascope_id_species,
                   file = file.path(out_dir,
                                    paste0(sample_name, ".metascope_species.csv")))
  message("Saving metascope grouped species to ",
          file.path(out_dir, paste0(sample_name, ".metascope_species.csv")))
  
  # Create fasta directory in tmp directory to save fasta sequences
  fastas_tmp_dir <- file.path(tmp_dir, "fastas")
  if(!dir.exists(fastas_tmp_dir)) dir.create(fastas_tmp_dir,
                                             recursive = TRUE)
  unlink(paste0(fastas_tmp_dir, "/*"), recursive = TRUE)
  
  message("Generating fasta sequences from bam file")
  # How many taxa
  num_taxa_loop <- min(nrow(metascope_id_species), num_results)
  # Extract sequence information from BAM file
  seq_info_df <- as.data.frame(Rsamtools::seqinfo(bam_file)) |>
    tibble::rownames_to_column("seqnames") |>
    dplyr::mutate("original_seqids" = !!dplyr::sym("seqnames"))
  
  # Obtain NCBI uids from readnames (accessions)
  if (db == "ncbi") {
    all_ids <- identify_rnames(list(list(rname = seq_info_df$seqnames))) |>
      taxize::genbank2uid(key = NCBI_key) |>
      lapply(function(x) x[[1]]) |> unlist()
    seq_info_df$seqnames <- all_ids
  }
  # Iterate over 
  write_fastas <- function(i) {
    taxids <- strsplit(metascope_id_species$TaxonomyIDs[i], split = ",")[[1]]
    read_proportions <- strsplit(metascope_id_species$read_proportions[i],
                                 split = ",")[[1]]
    reads_to_sample <- ceiling(as.numeric(read_proportions) * num_reads)
    if (length(reads_to_sample) > num_reads) {
      taxids <- taxids[seq_len(num_reads)]
      reads_to_sample <- reads_to_sample[seq_len(num_reads)]
    }
    ids_n <- lapply(seq_along(taxids), function(i) c(taxids[i], reads_to_sample[i]))
    # Parallelize for loop
    seqs_list <- BiocParallel::bplapply(ids_n, get_multi_seqs, bam_file = bam_file,
                                        seq_info_df = seq_info_df,
                                        metascope_id_tax = metascope_id_tax,
                                        sorted_bam_file = sorted_bam_file, BPPARAM = param)
    seqs <- do.call(c, seqs_list)
    Biostrings::writeXStringSet(
      seqs$seq, filepath = file.path(fastas_tmp_dir, paste0(sprintf("%05d", i), ".fa")))
  }
  
  # Parallelize the writing of fastas
  out <- BiocParallel::bplapply(seq_len(num_taxa_loop), write_fastas, BPPARAM = param)
  
  # Create blast directory in tmp directory to save blast results in
  blast_tmp_dir <- file.path(tmp_dir, "blast")
  if(!dir.exists(blast_tmp_dir)) dir.create(blast_tmp_dir, recursive = TRUE)
  unlink(paste0(blast_tmp_dir, "/*"), recursive = TRUE)
  
  # Create accessions database
  if (is.null(accessions_path)) {
    accessions_path <- file.path(tmp_dir, 'accessionTaxa.sql')
    taxonomizr::prepareDatabase(accessions_path)
  }
  
  # Run rBlast on all metascope microbes
  message("Running BLASTN on all sequences")
  blastn_results(results_table = metascope_id_species, bam_file = bam_file,
                 num_results = num_results, num_reads_per_result = num_reads,
                 hit_list = hit_list, num_threads = num_threads,
                 db_path = db_path, out_path = blast_tmp_dir,
                 sample_name = sample_name, quiet = TRUE,
                 accessions_path = accessions_path, fasta_dir = fastas_tmp_dir,
                 param = param)
  
  # Run Blast metrics
  message("Running BLAST metrics on all blast results")
  blast_result_metrics_df <- plyr::adply(
    list.files(blast_tmp_dir, full.names = TRUE), 1, blast_result_metrics,
    accessions_path = accessions_path, db = db)
  blast_result_metrics_df <- blast_result_metrics_df[ , -which(names(blast_result_metrics_df) %in% c("X1"))]
  
  # Append Blast Metrics to MetaScope results
  if (nrow(metascope_id_species) > nrow(blast_result_metrics_df)) {
    ind <- seq(nrow(blast_result_metrics_df) + 1, nrow(metascope_id_species))
    blast_result_metrics_df[ind, ] <- NA
  }
  print_file <- file.path(out_dir, paste0(sample_name, ".metascope_blast.csv"))
  
  metascope_blast_df <- data.frame(metascope_id_species, blast_result_metrics_df)
  utils::write.csv(metascope_blast_df, print_file)
  message("Results written to ", print_file)
  return(utils::head(metascope_blast_df))
}

#' Reassign reads from Blast alignment
#'
#' This function allows the user to reassign reads who's NCBI BLAST results
#' contradict results provided by MetaScope to assignments that were BLAST
#' validated. It returns an updated csv with reads reassigned according to their
#' BLAST validation.
#'
#' @param metascope_blast_df metascope blast dataframe; this is the output of metascope_blast
#' that was written to a CSV file.
#' @param species_threshold species threshold for validation
#' Species percentage hit (Blasted reads) threshold to validate a taxa
#' Threshold >0.5 to be considered  (50%, so this is a proportion, can be 0.2 default)
#' @param num_hits number of hits 
#' @param blast_tmp_dir directory with blast results
#' @inheritParams metascope_blast
#'
#' @returns Returns an updated blast reassignment file
#'
#' @export

blast_reassignment <- function(metascope_blast_df, species_threshold, num_hits,
                               blast_tmp_dir, out_dir, sample_name) {
  # Create validated column to determine if reads should be reassigned to accession
  metascope_blast_df <- metascope_blast_df |>
    dplyr::mutate(blast_validated = (!!dplyr::sym("species_percentage_hit") > species_threshold),
                  index = seq_len(nrow(metascope_blast_df)))
  blast_files <- list.files(blast_tmp_dir, full.names = TRUE)
  # If num_hits is small
  num_hits <- min(num_hits, nrow(metascope_blast_df))
  test_indices <- which(!metascope_blast_df$blast_validated)
  get_blast_summary <- function(i) {
    blast_summary <-
      tryCatch(
        {
          blast_summary <- readr::read_csv(blast_files[i],
                                           show_col_types = FALSE) |>
            dplyr::group_by(!!dplyr::sym("qseqid")) |>
            dplyr::slice_min(!!dplyr::sym("evalue"), with_ties = TRUE) |>
            # Removing duplicate query num and query species
            dplyr::distinct(!!dplyr::sym("qseqid"), !!dplyr::sym("species"), .keep_all = TRUE) |>
            dplyr::ungroup() |>
            dplyr::group_by( !!dplyr::sym("genus"), !!dplyr::sym("species")) |>
            dplyr::summarise("num_reads" = dplyr::n(), .groups="keep") |>
            dplyr::slice_max(order_by = !!dplyr::sym("num_reads"), with_ties = TRUE) |>
            dplyr::left_join(metascope_blast_df[1:num_hits, ],
                             by = dplyr::join_by(genus == best_hit_genus,
                                                 species == best_hit_species)) |>
            dplyr::filter(!!dplyr::sym("blast_validated") == TRUE) |>
            dplyr::mutate(reassignment_proportion = !!dplyr::sym("num_reads") / sum(!!dplyr::sym("num_reads")),
                          reassigned_read_count = metascope_blast_df$read_count[i] * !!dplyr::sym("reassignment_proportion"),
                          reassigned_Proportion = metascope_blast_df$Proportion[i] * !!dplyr::sym("reassignment_proportion"),
                          reassigned_readsEM = metascope_blast_df$readsEM[i] * !!dplyr::sym("reassignment_proportion"),
                          reassigned_EMProportion = metascope_blast_df$EMProportion[i] * !!dplyr::sym("reassignment_proportion"))
          blast_summary
        },
        error = function(cond) {
          message("Problem with blast summary for current file: ")
          message(blast_files[i])
          message(cond)
          # Choose a return value in case of error
          blast_summary <- dplyr::tibble(
            genus = numeric(),
            species = numeric(),
            num_reads = numeric()
          )
          blast_summary
        }
      )
    return(blast_summary)
  }
  
  all_blast_summaries <- plyr::alply(test_indices, 1,
                                     get_blast_summary) |>
    magrittr::set_names(paste0("X", test_indices))
  drop_indices_PRE <- plyr::laply(all_blast_summaries,
                                  function(x) nrow(x) > 0)
  drop_indices_PRE_2 <- test_indices[drop_indices_PRE]
  drop_indices <- drop_indices_PRE_2[drop_indices_PRE_2 != 1]
  
  all_names <- c("read_count", "Proportion", "readsEM", "EMProportion")
  summarize_unval <- function(i) {
    this_blast_summary <- all_blast_summaries[[paste0("X", i)]]
    
    summarize_reassignments <- function(n_ind) {
      metascope_index <- this_blast_summary$index[n_ind]
      init_list <- vector(mode = "list", length = length(all_names)) |>
        magrittr::set_names(all_names)
      output <- plyr::aaply(names(init_list), 1, function(x)
        metascope_blast_df[metascope_index, x] +
          this_blast_summary[n_ind, x]) |>
        magrittr::set_names(all_names) |> as.data.frame()
      return(output)
    }
    output <- plyr::adply(seq_len(nrow(this_blast_summary)), 1,
                          summarize_reassignments, .id = NULL) |>
      as.data.frame() |>
      dplyr::mutate(id = this_blast_summary$index)
    return(output)
  }
  to_replace <- plyr::adply(drop_indices, 1, summarize_unval, .id = NULL) |>
    dplyr::distinct() |> dplyr::arrange(id)
  reassigned_metascope_blast <- metascope_blast_df |> dplyr::select(-"index")
  reassigned_metascope_blast[to_replace$id, all_names] <- to_replace[, all_names] 
  if (length(drop_indices) > 0) {
    reassigned_metascope_blast <- reassigned_metascope_blast[-drop_indices, ]
  }
  # Clean up unused columns
  reassigned_metascope_blast <- reassigned_metascope_blast |>
    dplyr::select(-c("IDs", "TaxonomyIDs", "read_proportions",
                     "Proportion", "readsEM", "EMProportion", "X")) |>
    dplyr::rename('reassigned_read_count' = 'read_count')
  print_file <- file.path(out_dir, paste0(sample_name, ".metascope_blast_reassigned.csv"))
  utils::write.csv(reassigned_metascope_blast, print_file)
  message("Results written to ", print_file)
  return(reassigned_metascope_blast)
}

