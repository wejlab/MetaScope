getExtension <- function(file){
  ex <- strsplit(basename(file), split="\\.")[[1]]
  return(ex[-1])
}

# download_accessions.R
#' Download indexes required for MetaScope ID and MetaBlast
#'
#' This function will automatically download the NCBI accessions database, the
#' SILVA taxonomy database, and the NCBI Blast 16S database and prepare it for
#' MetaScope ID and MetaBlast
#'
#' @param ind_dir Character path to directory where indices should be saved
#' @param tmp_dir Character path to directory for storing temp files. (Useful
#' to avoid redownloading) Defaults to \code{file_path(ind_dir, "tmp")}
#' @param remove_tmp_dir Delete tmp_dir after downloads are complete? Defaults
#' to \code{TRUE}
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
#' @export

download_accessions <- function(ind_dir,
                                tmp_dir = file_path(ind_dir, "tmp"),
                                remove_tmp_dir = TRUE,
                                NCBI_accessions_database = TRUE,
                                NCBI_accessions_name = "accessionTaxa",
                                silva_taxonomy_database = TRUE,
                                silva_taxonomy_name = "all_silva_headers",
                                blast_16S_database = TRUE,
                                blast_16S_name = "16S_ribosomal_RNA") {
  # Check if ind_dir exists
  if (!file.exists(ind_dir)) {
    stop("ind_dir is not a valid directory")
  }
  if(!file.exists(tmp_dir)){
    dir.create(tmp_dir)
  }
  if (NCBI_accessions_database) {
    message("Downloading NCBI accessions database")
    if (!identical(getExtension(NCBI_accessions_name), "sql")) {
      NCBI_accessions_name <- paste0(NCBI_accessions_name, ".sql")
    }
    taxonomizr::prepareDatabase(file.path(ind_dir, NCBI_accessions_name),
                                tmpDir = tmp_dir)
  }
  if (silva_taxonomy_database) {
    message("Downloading silva taxonomy database")
    location <- "https://github.com/wejlab/metascope-docs/raw/main/all_silva_headers.rds"
    if (!identical(getExtension(silva_taxonomy_name), "rds")) {
      NCBI_accessions_name <- paste0(silva_taxonomy_name, ".rds")
    }
    destination <- paste(ind_dir, silva_taxonomy_name, sep = "/")
    utils::download.file(location, destination)
  }
  if (blast_16S_database) {
    message("Downloading Blast 16S rRNA database")
    location <- "https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz"
    destination <- paste(ind_dir, blast_16S_name, sep = "/")
    utils::download.file(location, destination)
  }
  if (remove_tmp_dir) {
    unlink(tmp_dir)
  }
}

