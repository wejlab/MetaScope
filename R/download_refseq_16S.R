#' Download RefSeq 16S rRNA Bacterial and Archael libraries
#'
#' This function will automatically download the 16S rRNA RefSeq libraries from 
#' \code{https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/} and combine them 
#' into a single .fna file
#'
#' @param out_dir Character string giving the name of the directory to which
#'   libraries should be output. Defaults to creation of a new temporary
#'   directory.
#' @param combined_name Name of output combined file. Default is 
#'   \code{"refseq_16S.fna"}
#'
#' @return Returns a character string of the path to the combined 16S rRNA .fna 
#'    file
#'    
#' @note This function requires the suggested packages \pkg{RCurl} and 
#'    \pkg{R.utils}. If they are not installed, you will need to install 
#'    them manually.
#'
#' @export
#' @examples
#' #### Download 16S rRNA Genomes
#' 
#' download_refseq_16S <- function(out_dir = "out_dir", 
#'   combined_name = "refseq_16S.fna")
#' 


download_refseq_16S <- function(out_dir,
                                combined_name = "refseq_16S.fna") {
  # Check for required suggested packages
  if (!requireNamespace("RCurl", quietly = TRUE)) {
    stop("Package 'RCurl' is required for download_refseq_16s_only(). Please install it.", call. = FALSE)
  }
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    stop("Package 'R.utils' is required for download_refseq_16s_only(). Please install it.", call. = FALSE)
  }
  
  # Import namespaces explicitly
  RCurl <- asNamespace("RCurl")
  R.utils <- asNamespace("R.utils")
  
  # Create Directory
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  base_urls <- c(
    "ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/",
    "ftp://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/"
  )
  
  download_from_url <- function(url, local_dir) {
    listing <- str_split(getURL(url, ftp.use.epsv = TRUE, dirlistonly = TRUE), "\n")[[1]]
    listing <- str_trim(listing)
    listing <- listing[listing != ""]
    
    for (item in listing) {
      if (grepl("/$", item)) {
        # Recurse into subdirectories
        new_url <- paste0(url, item)
        download_from_url(new_url, local_dir)
      } else if (grepl("16S.*\\.fna\\.gz$", item)) {
        # Only download 16S FASTA files
        file_url <- paste0(url, item)
        dest_file <- file.path(local_dir, item)
        tryCatch({
          download.file(file_url, destfile = dest_file, quiet = TRUE, mode = "wb")
          message("Downloaded: ", file_url)
        }, error = function(e) {
          message("Failed to download: ", file_url)
        })
      }
    }
  }
  
  # Download all 16S files
  for (base_url in base_urls) {
    download_from_url(base_url, out_dir)
  }
  
  # Gunzip all downloaded files
  gz_files <- list.files(out_dir, pattern = "\\.fna\\.gz$", full.names = TRUE)
  message("Decompressing ", length(gz_files), " files...")
  
  for (gz in gz_files) {
    tryCatch({
      gunzip(gz, overwrite = TRUE, remove = TRUE)  # decompress & remove .gz
      message("Unzipped: ", gz)
    }, error = function(e) {
      message("Failed to unzip: ", gz)
    })
  }
  
  # Concatenate into one combined FASTA
  fna_files <- list.files(out_dir, pattern = "\\.fna$", full.names = TRUE)
  combined_path <- file.path(out_dir, combined_name)
  
  message("Concatenating into: ", combined_path)
  file.create(combined_path)  # ensure file exists
  for (fna in fna_files) {
    system(paste("cat", shQuote(fna), ">>", shQuote(combined_path)))
  }
  
  message("All 16S FASTAs combined into: ", combined_path)
  
  # Delete original fna_files
  unlink(fna_files)
  
  return(combined_path)
}