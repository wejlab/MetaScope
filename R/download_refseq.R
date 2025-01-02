# Helper function to identify parent kingdom, rank
id_kingdom_rank <- function(classification_table, taxon, rank_input, quiet) {
  # Get the parent taxon in super kingdom rank
  parent_kingdom <- classification_table$name[
    classification_table$rank == "superkingdom"]
  parent_rank <- "Superkingdom"
  # The Eukaryota superkingdom isn't a folder on the NCBI FTP site
  # Instead, grab the corresponding kingdom
  # because some eukaryota kingdoms are folders on the FTP site
  if (identical(parent_kingdom, "Eukaryota")) {
    parent_kingdom <- classification_table$name[
      classification_table$rank == "kingdom"]
    parent_rank <- "Kingdom"
    # There are only two eukaryota kingdoms in the taxonomy table that
    # have a corresponding folder on the FTP site. If not one of these
    # two kingdoms then forced to search all the eukaryota folders on the
    # site
    if (length(parent_kingdom) != 0) {
      if (!(parent_kingdom %in% c("Viridiplantae", "Fungi"))) {
        parent_kingdom <- classification_table$name[
          classification_table$rank == "superkingdom"]
        parent_rank <- "Superkingdom"
      }
    } else {
      parent_kingdom <- classification_table$name[
        classification_table$rank == "superkingdom"]
      parent_rank <- "Superkingdom"
    }
  }
  if (!quiet) message(taxon, " is a ", rank_input, " under the ",
                      parent_kingdom, " ", parent_rank)
  parent_kingdom <- tolower(parent_kingdom)
  # Some folders on the NCBI FTP site don't match the parent kingdom taken
  # from the classification table. We'll apply some modifications.
  if (identical(parent_kingdom, "viruses")) parent_kingdom <- "viral"
  if (identical(parent_kingdom, "viridiplantae")) parent_kingdom <- "plant"
  return(parent_kingdom)
}

# Helper function to download the refseq table (for a given group)
get_table <- function(what_tax) {
  refseq_link <- paste("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/",
                       what_tax, "/assembly_summary.txt", sep = "")
  refseq_table <- utils::read.table(refseq_link, header = TRUE,
                                    sep = "\t", comment.char = "",
                                    quote = "", skip = 1, fill = TRUE)
  return(refseq_table)
}

# Helper function to obtain refseq table
download_parentkingdom <- function(parent_kingdom, quiet) {
  # Download assembly summary refseq table from NCBI
  # (includes genome download link)
  if (!quiet) message("Loading the refseq table for ", parent_kingdom)
  to_grab <- c("archaea", "bacteria", "fungi", "plant", "viral")
  if (parent_kingdom %in% to_grab) {
    refseq_table <- get_table(parent_kingdom)
    return(refseq_table)
  } else if (parent_kingdom == "eukaryota") {
    refseq_table <- lapply(c("vertebrate_mammalian", "vertebrate_other",
                             "invertebrate", "protozoa"), get_table) |>
      data.table::rbindlist() |> as.data.frame()
    return(refseq_table)
  } else message("Parent kingdom table could not be retrieved from",
                 "NCBI database. Try a different taxon.")
}

# Helper function to get species table
get_speciestab <- function(children_list, refseq_table, taxon,
                           representative, reference, quiet) {
  if (!quiet) message("Creating table of relevant taxa")
  # Filter table, keep lines with species or strains of input
  # If the taxon specified is a strain or species with no strain,
  # filter table by the taxon
  # Otherwise, filter table by the taxon in the children list
  if (rlang::is_empty(children_list)) {
    species_table <- refseq_table[which(tolower(
      refseq_table$organism_name) %in% tolower(taxon)), ]
  } else {
    species_table <- refseq_table[which(tolower(
    refseq_table$organism_name) %in% tolower(children_list)), ]
  }
  # Reduce table size based on reference or representative
  if (representative) reference <- TRUE
  if (representative && reference) {
    species_table <- species_table[
      species_table$refseq_category %in% c("reference genome",
                                           "representative genome"), ]
  } else if (!representative && reference) {
    species_table <- species_table[
      species_table$refseq_category == "reference genome", ]
  }
  return(species_table)
}

# Helper function to get cache for BiocFileCache
.get_cache <- function() {
  cache <- tools::R_user_dir("MetaScope", which = "cache")
  BiocFileCache::BiocFileCache(cache)
}


# Another helper function to download genomes
download_genomes <- function(species_table, taxon, patho_out, compress,
                             out_dir, caching, quiet) {
  total_genomes <- nrow(species_table)
  if (total_genomes == 0) {
    stop("No available genome for ", taxon,
         " - try setting both `representative` and `reference`",
         " to FALSE to obtain genome from the nucleotide database.")
  } else {
    if (!quiet) message("Downloading ", total_genomes, " ", taxon,
                        " genome(s) from NCBI")
    # Delete existing genome files and combined fasta
    taxon <- taxon %>%
      stringr::str_replace_all("/", "_") %>%
      stringr::str_replace_all(" ", "_")
    download_dir <- file.path(out_dir, paste(taxon, "refseq_download",
                                             sep = "_"))
    # Remove any existing files/directories
    if (file.exists(download_dir)) {
      unlink(file.path(download_dir, "*")) %>% message()
      unlink(download_dir, recursive = TRUE, force = TRUE) %>% message()
    }
    if (compress) {
      combined_fasta <- file.path(out_dir, paste(taxon, "fasta.gz", sep = "."))
      combined_fasta_patho <- file.path(out_dir,
                                        paste(taxon, "pathoscope.fasta.gz",
                                              sep = "."))
    } else {
      combined_fasta <- file.path(out_dir, paste(taxon, "fasta", sep = "."))
      combined_fasta_patho <- file.path(
        out_dir, paste(taxon, "pathoscope.fasta", sep = "."))
    }
    # Start with a new combined file
    if (file.exists(combined_fasta)) unlink(combined_fasta, force = TRUE)
    if (file.exists(combined_fasta_patho)) unlink(combined_fasta_patho,
                                                  force = TRUE)
  }
  # Download the genome
  for (i in seq_len(nrow(species_table))) {
    tryCatch({
      genome_file <- paste(basename(as.character(
        species_table[i, ]$ftp_path)), "genomic.fna.gz", sep = "_")
      location <- paste(species_table[i, ]$ftp_path, genome_file,
                        sep = "/")
      if (!caching) {
        if (!dir.exists(download_dir)) dir.create(download_dir)
        destination <- paste(download_dir, genome_file, sep = "/")
        utils::download.file(location, destination)
      } else if (caching) {
        bfc <- .get_cache()
        rid <- BiocFileCache::bfcquery(bfc, genome_file, "rname")$rid
        if (!length(rid)) {
          if (i %% 10 == 0 && !quiet) {
            message("Number of Genomes Downloaded: ", i, "/",
                    total_genomes, " (",
                    round(100 * i / total_genomes, 2), "%)")
            }
          rid <- names(BiocFileCache::bfcadd(bfc, genome_file, location))
        }
        if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid))) {
          BiocFileCache::bfcdownload(bfc, rid)
          BiocFileCache::bfcrpath(bfc, rids = rid)
        } else {
          if (!quiet) message("Caching is set to TRUE, ",
                              "and it appears that this file is already downloaded ",
                              "in the cache. It will not be downloaded again.")
        }
        destination <- BiocFileCache::bfcrpath(bfc, rids = rid)
      }
      # Write genome to concatenated file
      in_con <- file(destination, open = "rb")
      ref <- gzcon(in_con) %>% base::readLines()
      close(in_con)
      ref %>% data.table::as.data.table() %>%
        data.table::fwrite(file = combined_fasta, append = TRUE,
                           quote = FALSE, sep = " ", compress = "gzip",
                           col.names = FALSE, row.names = FALSE)
      # Format for pathoscope and write to file
      if (patho_out) {
        # Identify accessions
        ind <- ref %>% stringr::str_starts(">")
        accession <- ref[ind] %>% stringr::str_split(pattern = " ") %>%
          vapply(S4Vectors::head, n = 1, FUN.VALUE = character(1)) %>%
          stringr::str_remove_all(">")
        ref[ind] <- paste(">ti|", species_table[i, ]$taxid, "|org|",
                          gsub(" ", "_", species_table[i, ]$organism_name),
                          "|accession|", accession, sep = "")
        ref %>% data.table::as.data.table() %>%
          data.table::fwrite(file = combined_fasta_patho,
                             append = TRUE, quote = FALSE, sep = " ",
                             compress = "gzip", col.names = FALSE,
                             row.names = FALSE)
      }
    }, error = function(e) message(conditionMessage(e)))
  }
  unlink(file.path(download_dir, "*"), force = TRUE)
  unlink(file.path(download_dir), recursive = TRUE, force = TRUE)
  if (!quiet) message("DONE! ", i, " genome(s) saved to file ",
                      file.path(combined_fasta))
  return(combined_fasta)
}

find_strains <- function(intable) {
  all_rank <- which(intable$rank == "no rank")
  if (all(is.na(all_rank))) return(intable)
  for (t_n in all_rank[all_rank > 1]) {
    if (intable$rank[t_n] == "no rank" && intable$rank[t_n - 1] == "species") {
      intable$rank[t_n] <- "strain"
    }
  }
  return(intable)
}

#' Download RefSeq genome libraries
#'
#' This function will automatically download RefSeq genome libraries in a fasta
#' format from the specified taxon. The function will first download the
#' summary report at:
#' \code{ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/**kingdom**/assembly_summary.txt},
#' and then use this file to download the genome(s) and combine them in a
#' single compressed or uncompressed .fasta file.
#'
#' When selecting the \code{taxon} to be downloaded, if you receive an error
#' saying \code{Your input is not a valid taxon}, please take a look at the
#' \code{taxonomy_table} object, which can be accessed with the command
#' \code{MetaScope:::taxonomy_table)}. Only taxa with exact spelling as they appear
#' at any level of the table will be acknowledged.
#'
#' @param taxon Name of single taxon to download. The taxon name should be a
#'   recognized NCBI scientific or common name, with no grammatical or
#'   capitalization inconsistencies. All available taxonomies are visible by
#'   accessing the \code{MetaScope:::taxonomy_table} object included in the package.
#' @param reference Download only RefSeq reference genomes? Defaults to
#'   \code{TRUE}. Automatically set to \code{TRUE} if \code{representative =
#'   TRUE}.
#' @param representative Download RefSeq representative and reference genomes?
#' Defaults to \code{FALSE}. If \code{TRUE}, reference is automatically set at
#'   \code{TRUE}.
#' @param compress Compress the output .fasta file? Defaults to \code{TRUE}.
#' @param patho_out Create duplicate outpute files compatible with PathoScope?
#'   Defaults to \code{FALSE}.
#' @param out_dir Character string giving the name of the directory to which
#'   libraries should be output. Defaults to creation of a new temporary
#'   directory.
#' @param caching Whether to use BiocFileCache when downloading genomes.
#'   Default is \code{FALSE}.
#' @param accessions_path (character) Filepath to NCBI accessions SQL
#'   database. See \code{taxonomzr::prepareDatabase()}.
#' @param quiet Turns off most messages. Default is \code{TRUE}.
#'
#' @return Returns a .fasta or .fasta.gz file of the desired RefSeq genomes.
#' This file is named after the kingdom selected and saved to the current
#' directory (e.g. 'bacteria.fasta.gz'). This function also has the option
#' to return a .fasta file formatted for PathoScope as well
#' (e.g. bacteria.pathoscope.fasta.gz') if \code{path_out = TRUE}.
#'
#' @export
#' @examples
#' #### Download RefSeq genomes
#'
#' ## Download all RefSeq reference Bovismacovirus genus genomes
#' download_refseq('Bovismacovirus', reference = FALSE, representative = FALSE,
#'                 out_dir = NULL, compress = TRUE, patho_out = FALSE,
#'                 caching = TRUE)
#'

download_refseq <- function(taxon, reference = TRUE, representative = FALSE,
                            compress = TRUE, patho_out = FALSE,
                            out_dir = NULL, caching = FALSE,
                            quiet = TRUE, accession_path = NULL) {
  if (is.null(out_dir)) {
    out_dir <- tempfile()
    dir.create(out_dir)
  }
  if (is.null(accession_path)) {
    stop("Process halted. Requires an taxonomizr accession_path")
  }
  taxon_id <- taxonomizr::getId(taxon, sqlFile = accession_path)
  classification_table <- taxonomizr::getRawTaxonomy(taxon_id, accession_path)
  classification_table <- classification_table[[1]] |> as.data.frame() |>
    tibble::rownames_to_column(var = "rank")
  colnames(classification_table) <- c("rank", "name")
  rank_input <- classification_table$rank[1]
  parent_kingdom <- id_kingdom_rank(classification_table, taxon, rank_input, quiet)
  refseq_table <- download_parentkingdom(parent_kingdom, quiet)
  # Get NCBI scientific names of children species or strains
  taxonomy_table <- get0("taxonomy_table", envir = asNamespace("MetaScope"))
  if (rank_input == "no rank") base::stop("No rank detected")
  children_list <- get_children(taxon, rank_input, tax_dat = taxonomy_table)
  species_table <- get_speciestab(children_list, refseq_table, taxon,
                                  representative, reference, quiet)
  combined_fasta <- download_genomes(species_table, taxon, patho_out,
                                     compress, out_dir, caching, quiet)
  return(combined_fasta)
}