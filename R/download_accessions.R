# download_accessions.R
#' Download indexes required for MetaScope ID and MetaBlast
#'
#' This function will automatically download the NCBI accessions database, the
#' SILVA taxonomy database, and the NCBI Blast 16S database and prepare it for
#' MetaScope ID and MetaBlast
#'
#' @param ind_dir Character path to directory where indices should be saved
#' @param NCBI_accessions_database Download taxonomizr NCBI accessions database?
#' Defaults to \code{TRUE}.
#' @param NCBI_accessions_name Name to save taxonomizr NCBI accessions database
#' @param silva_taxonomy_database Download silva taxonomy database? Defaults to
#' \code{TRUE}.
#' @param silva_taxonomy_name  Name of silva headers file. Defaults
#' to the one hosted in MetaScope, \code{"all_silva_headers.rds"}.
#' @param NCBI_accessions_name Name to save silva taxonomy database
#' @param blast_16S_database Download NCBI 16S Blast database? Defaults to
#' \code{TRUE}.
#' @param blast_16S_name Name to save NCBI 16S Blast database

download_accessions <- function(ind_dir,
                             NCBI_accessions_database = TRUE,
                             NCBI_accessions_name = "accessionTaxa.sql",
                             silva_taxonomy_database = TRUE,
                             silva_taxonomy_name = "all_silva_headers.rds",
                             blast_16S_database = TRUE,
                             blast_16S_name = "16S_ribosomal_RNA") {
  if (NCBI_accessions_database) {
    message("Downloading NCBI accessions database")
    taxonomizr::prepareDatabase(file.path(ind_dir, NCBI_accessions_name))
  }
  if (silva_taxonomy_database) {
    message("Downloading silva taxonomy database")
    location <- "https://github.com/wejlab/metascope-docs/raw/main/all_silva_headers.rds"
    destination <- paste(ind_dir, silva_taxonomy_name, sep = "/")
    utils::download.file(location, destination)
  }
  if (blast_16S_database) {
    message("Downloading Blast 16S rRNA database")
    location <- "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz"
    destination <- paste(ind_dir, blast_16S_name, sep = "/")
    utils::download.file(location, destination)
  }
}

