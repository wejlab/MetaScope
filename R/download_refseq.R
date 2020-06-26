#' Download RefSeq genome libraries
#'
#' This function will automatically download RefSeq genome libraries in a
#' .fasta format from the specified kingdom. The function will first
#' download the summary report at:
#' ftp://ftp.ncbi.nlm.nih.gov/genomes/**kingdom**/overview.txt,
#' and then use this file to download genomes and combine them in a single
#' compressed or uncompressed .fasta file.
#' @param kingdom Select the kingdom taxonomy to download.
#' Options are \code{'archaea'}, \code{'bacteria'}, \code{'fungi'},
#' \code{'invertebrate'}, \code{'plant'}, \code{'protozoa'},
#' \code{'vertebrate_mammalian'}, \code{'vertebrate_other'}, or \code{'viral'}.
#' @param reference Download only RefSeq reference genomes?
#' Defaults to \code{TRUE}.
#' Automatically set to \code{TRUE} if \code{representative = TRUE}.
#' @param representative Download only RefSeq representative genomes?
#' Defaults to \code{FALSE}.
#' If \code{TRUE}, reference is automatically set at \code{TRUE}.
#' @param compress Compress the output .fasta file? Defaults to \code{TRUE}.
#' @param patho_out Create duplicate outpute files compatible with PathoScope?
#' Defaults to \code{FALSE}.
#' @return Returns a .fasta or .fasta.gz file of the desired RefSeq genomes.
#' This file is named after the kindom selectd and saved to the current
#' directory (e.g. 'bacteria.fasta.gz'). Currently, this function also returns
#' a .fasta file formatted for PathoScope as well
#' (e.g. bacteria.pathoscope.fasta.gz') if \code{path_out = TRUE}.
#'
#' @examples
#' ## Download all RefSeq reference bacterial genomes
#' download_refseq('bacteria')
#'
#' ## Download all RefSeq representative viral genomes
#' download_refseq( 'viral', representative = TRUE )
#'
#' ## Download all RefSeq viral genomes
#' download_refseq( 'viral', reference = FALSE )
#'
#' @export
#' 

download_refseq <- function(kingdom, reference = TRUE,
                            representative = FALSE,
                            compress = TRUE, patho_out = FALSE) {
  ## check if user provided a valid kingdom
  kingdom_list <- c("archaea", "bacteria", "fungi", "invertebrate", "plant",
                    "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral")
  if (!(kingdom %in% kingdom_list)) {
    stop("You supplied a kingdom not in the kingdom list")
  }

  ## Download refseq table
  table_name <- paste("ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/", kingdom,
                      "/assembly_summary.txt", sep = "")
  message(paste("Downloading kingdom assembly list from:", table_name))
  kingdom_table <- read.table(table_name, header = T, sep = "\t",
                              comment.char = "",
                              quote = "", skip = 1)

  ## Reduce the table size based on reference or represenative
  if (representative) {
    reference <- TRUE
  }
  if (representative & reference) {
    king_table <- kingdom_table[kingdom_table$refseq_category %in%
                                  c("reference genome",
                                    "representative genome"), ]
  } else {
    if (!representative & reference) {
      king_table <- kingdom_table[kingdom_table$refseq_category ==
                                    "reference genome", ]
    } else {
      king_table <- kingdom_table
    }
  }
  total_genomes <- nrow(king_table)
  message(paste("Downloading", total_genomes, kingdom, "genomes from RefSeq"))

  ## delete existing genome files and combined fasta--make these user-defined

  download_dir <- paste(kingdom, "refseq_download", sep = "_")
  # remove any existing files/directories
  unlink(download_dir, recursive = TRUE, force = TRUE)
  if (compress) {
    combined_fasta <- paste(kingdom, "fasta.gz", sep = ".")
    combined_fasta_patho <- paste(kingdom, "pathoscope.fasta.gz", sep = ".")
  } else {
    combined_fasta <- paste(kingdom, "fasta", sep = ".")
    combined_fasta_patho <- paste(kingdom, "pathoscope.fasta", sep = ".")
  }
  tryCatch({
    suppressWarnings(file.remove(combined_fasta))
    suppressWarnings(file.remove(combined_fasta_patho))
  })  # start with a new combined file

  ## Download genomes
  for (i in 1:nrow(king_table)) {
    tryCatch({
      if (i%%10 == 0) {
        message(paste("Number of Genomes Downloaded: ", i, "/",
                      total_genomes, " (", round(100 * i/total_genomes, 2),
                      "%)", sep = ""))
      }

      ## Download the genome
      genome_file <- paste(basename(as.character(king_table[i, ]$ftp_path)),
                           "genomic.fna.gz", sep = "_")
      location <- paste(king_table[i, ]$ftp_path, genome_file, sep = "/")
      destination <- paste(download_dir, genome_file, sep = "/")
      if (!dir.exists(download_dir)) {
        dir.create(download_dir)
      }
      download.file(location, destination)

      ## read in the genome
      ref <- Biostrings::readDNAStringSet(destination)

      ## write to file
      Biostrings::writeXStringSet(ref, combined_fasta, append = T,
                                  compress = compress)

      ## format for pathoscope and write to file
      accession <- NULL
      for (j in strsplit(names(ref), " ")) {
        accession <- c(accession, j[1])
      }
      names(ref) <- paste("ti|", king_table[i, ]$taxid, "|org|",
                          gsub(" ", "_", king_table[i, ]$organism_name),
                          "|accession|", accession, sep = "")
      Biostrings::writeXStringSet(ref, combined_fasta_patho, append = T,
                                  compress = compress)

      ## delete intermediate download files
      unlink(download_dir, recursive = TRUE, force = TRUE)
    }, error = function(e) {
      cat("ERROR :", conditionMessage(e), "\n")
    })
  }
  # Remove pathoscope file if unwanted
  if (!patho_out) file.remove(combined_fasta_patho)

  # Ensure removal of intermediate folder of files
  unlink("viral_refseq_download", recursive = TRUE, force = TRUE)

  message("DONE! Downloaded", i, "genomes to", combined_fasta, sep = " ")
}
