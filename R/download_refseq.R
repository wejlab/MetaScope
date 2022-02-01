globalVariables(c("taxonomy_table"))

#' Download RefSeq genome libraries
#'
#' This function will automatically download RefSeq genome libraries in a
#' .fasta format from the specified taxon. The function will first
#' download the summary report at:
#' ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/**kingdom**/assembly_summary.txt,
#' and then use this file to download the genome(s) and combine them in a single
#' compressed or uncompressed .fasta file.
#'
#' @param taxon Select one taxon to download. The taxon name should be a
#' recognized NCBI scientific name, with no grammatical or capitalization
#' inconsistencies. All available taxonomies are visible by accessing the
#' \code{taxonomy_table} object included in the package.
#' @param reference Download only RefSeq reference genomes?
#' Defaults to \code{TRUE}.
#' Automatically set to \code{TRUE} if \code{representative = TRUE}.
#' @param representative Download only RefSeq representative genomes?
#' Defaults to \code{FALSE}. If \code{TRUE}, reference is automatically
#' set at \code{TRUE}.
#' @param compress Compress the output .fasta file? Defaults to \code{TRUE}.
#' @param patho_out Create duplicate outpute files compatible with PathoScope?
#' Defaults to \code{FALSE}.
#'
#' @return Returns a .fasta or .fasta.gz file of the desired RefSeq genomes.
#' This file is named after the kindom selectd and saved to the current
#' directory (e.g. 'bacteria.fasta.gz'). Currently, this function also returns
#' a .fasta file formatted for PathoScope as well
#' (e.g. bacteria.pathoscope.fasta.gz') if \code{path_out = TRUE}.
#'
#' @export
#'
#' @examples
#' #### Download RefSeq genomes
#' 
#' ## Download all RefSeq reference bacterial superkingdom genomes
#' download_refseq('bacteria', reference = TRUE, representative = FALSE)
#'
#' ## Download all RefSeq representative mononegavirales genomes
#' download_refseq('mononegavirales', representative = TRUE)
#'
#' ## Download all RefSeq morbillivirus genomes
#' download_refseq('morbillivirus', reference = FALSE)
#'
#' ## Download all RefSeq bacilli reference genomes, uncompressed
#' download_refseq('Bacilli', reference = TRUE,
#'                 representative = FALSE, compress = FALSE)
#'
#' ## Download RefSeq Escherichia coli IAI1 strain 
#' download_refseq('Escherichia coli IAI1', reference = FALSE, compress = FALSE)
#' 

download_refseq <- function(taxon, reference = TRUE, representative = FALSE,
                            compress = TRUE, patho_out = FALSE) {
    
    # Get the rank of the input taxon
    tryCatch({suppressMessages(classification.table <- taxize::classification(
        taxize::get_uid(taxon, messages = FALSE)[[1]], db = 'ncbi')[[1]])},
        warning = function(w) {stop("Your input is not a valid taxon")}
    )
    message("Finding ", taxon)
    rank_input <- classification.table$rank[nrow(classification.table)]
    
    # Get the NCBI scientific names of children species or strains
    children_list <- get_children(taxon, rank_input, data = taxonomy_table)
    # Get the parent taxon in superkingdom rank
    parent_kingdom <- classification.table$name[classification.table$rank
                                                == "superkingdom"]
    parent_rank <- "Superkingdom"
    
    # The Eukaryota superkingdom is not a folder on the NCBI FTP site.
    # Instead grab the corresponding kingdom as some eukaryota kingdoms
    # are folders on the FTP site
    
    if (identical(parent_kingdom, "Eukaryota")){
        parent_kingdom <- classification.table$name[classification.table$rank
                                                    == "kingdom"]
        parent_rank <- "Kingdom"
        
        # There are only two eukaryota kingdoms in the taxonomy table that 
        # have a corresponding folder on the FTP site. If not one of these
        # two kingdoms then forced to search all the eukaryota folders on the
        # site
        if (length(parent_kingdom) != 0){
            if (!(parent_kingdom %in% c("Viridiplantae", "Fungi"))){
                parent_kingdom <- classification.table$name[classification.table$rank
                                                            == "superkingdom"]
                parent_rank <- "Superkingdom"
            }
        }
        else{
            parent_kingdom <- classification.table$name[classification.table$rank
                                                        == "superkingdom"]
            parent_rank <- "Superkingdom"
        }
    }
    
    message(taxon," is a ", rank_input, " under the ", parent_kingdom, " ",
            parent_rank)
    parent_kingdom <- tolower(parent_kingdom)
    
    # Some folders on the NCBI FTP site do not match the parent kingdom taken
    # from the classification table. Have to apply some modificiations.
    
    if (identical(parent_kingdom, "viruses")) parent_kingdom <- "viral"
    if (identical(parent_kingdom, "viridiplantae")) parent_kingdom <- "plant"

    # Download assembly summary refseq table from NCBI
        # (includes genome download link)
    message("Loading the refseq table for ", parent_kingdom)
    to_grab <- c("archaea", "bacteria", "fungi", "plant", "viral")
    if (parent_kingdom %in% to_grab) {
        refseq_link <- paste("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/",
                             parent_kingdom, "/assembly_summary.txt", sep = "")
        refseq_table <- utils::read.table(refseq_link, header = TRUE,
                                          sep = "\t", comment.char = "",
                                          quote = "", skip = 1)
    } else if (parent_kingdom == "eukaryota") {
        refseq_link_1 <- paste("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/",
                               "vertebrate_mammalian", "/assembly_summary.txt",
                               sep = "")
        refseq_table_1 <- utils::read.table(refseq_link_1, header = TRUE,
                                            sep = "\t", comment.char = "",
                                            quote = "", skip = 1)
        refseq_link_2 <- paste("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/",
                               "vertebrate_other", "/assembly_summary.txt",
                               sep = "")
        refseq_table_2 <- utils::read.table(refseq_link_2, header = TRUE,
                                            sep = "\t", comment.char = "",
                                            quote = "", skip = 1)
        refseq_link_3 <- paste("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/",
                               "invertebrate", "/assembly_summary.txt",
                               sep = "")
        refseq_table_3 <- utils::read.table(refseq_link_3, header = TRUE,
                                            sep = "\t", comment.char = "",
                                            quote = "", skip = 1)
        refseq_link_4 <- paste("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/",
                               "protozoa", "/assembly_summary.txt", sep = "")
        refseq_table_4 <- utils::read.table(refseq_link_4, header = TRUE,
                                            sep = "\t", comment.char = "",
                                            quote = "", skip = 1)
        refseq_table <- rbind(refseq_table_1, refseq_table_2,
                              refseq_table_3, refseq_table_4)
    } else message("Parent kingdom table could not be retrieved from",
                   "NCBI database. Try a different taxon.")
    
    # Filter the table, keep the lines with species or strains of input
    
    # If the taxon specifiec is a strain or species with no strain,
        # filter the table by the taxon
    # Otherwise filter the table by the taxon in the children list
    if (rlang::is_empty(children_list)) {
        species_table <- refseq_table[which(tolower(
            refseq_table$organism_name) %in% tolower(taxon)), ]
    } else species_table <- refseq_table[which(tolower(
        refseq_table$organism_name) %in% tolower(children_list)), ]

    # Reduce the table size based on reference or representative
    if (representative) reference <- TRUE
    if (representative & reference) {
        species_table <- species_table[
            species_table$refseq_category %in% c("reference genome",
                                                 "representative genome"), ]
    } else if (!representative & reference) species_table <- 
        species_table[species_table$refseq_category == "reference genome", ]

    total_genomes <- nrow(species_table)
    if (total_genomes == 0) {
        stop("No available genome for ", taxon, 
        " - try setting both `representative` and `reference`",
        " to either TRUE or FALSE")
    } else {
        message("Downloading ", total_genomes,
                " ", taxon, " genome(s) from RefSeq")
        # Delete existing genome files and combined fasta
        download_dir <- paste(taxon, "refseq_download", sep = "_")
        # Remove any existing files/directories
        unlink(download_dir, recursive = TRUE, force = TRUE)
        if (compress) {
            combined_fasta <- paste(taxon, "fasta.gz", sep = ".")
            combined_fasta_patho <- paste(taxon, "pathoscope.fasta.gz",
                                          sep = ".")
        } else {
            combined_fasta <- paste(taxon, "fasta", sep = ".")
            combined_fasta_patho <- paste(taxon, "pathoscope.fasta",
                                          sep = ".")
        }
        # Start with a new combined file
        tryCatch({
            suppressWarnings(file.remove(combined_fasta))
            suppressWarnings(file.remove(combined_fasta_patho))
        })
        ## Download the genome
        for (i in seq_len(nrow(species_table))) {
            tryCatch({
                if (i%%10 == 0) {
                    message("Number of Genomes Downloaded: ", i, "/",
                            total_genomes, " (",
                            round(100 * i/total_genomes, 2), "%)")
                }
                ## Download the genome
                genome_file <- paste(basename(as.character(
                    species_table[i, ]$ftp_path)), "genomic.fna.gz", sep = "_")
                location <- paste(species_table[i, ]$ftp_path, genome_file,
                                  sep = "/")
                destination <- paste(download_dir, genome_file, sep = "/")
                if (!dir.exists(download_dir)) dir.create(download_dir)
                utils::download.file(location, destination)
                ## Read in the genome
                ref <- Biostrings::readDNAStringSet(destination)
                ## Write to file
                Biostrings::writeXStringSet(ref, combined_fasta, append = TRUE,
                                            compress = compress)
                ## Format for pathoscope and write to file
                if (patho_out) {
                    accession <- NULL
                    for (j in strsplit(names(ref), " ")) {
                        accession <- c(accession, j[1])
                    }
                    names(ref) <- paste(
                        "ti|", species_table[i, ]$taxid, "|org|",
                        gsub(" ", "_", species_table[i, ]$organism_name),
                        "|accession|", accession, sep = "")
                    Biostrings::writeXStringSet(ref, combined_fasta_patho,
                                                append = TRUE,
                                                compress = compress)
                }
                ## Delete intermediate download files
                unlink(download_dir, recursive = TRUE, force = TRUE)
                }, error = function(e) cat("ERROR :", 
                                           conditionMessage(e), "\n"))
        }
        # Ensure removal of intermediate folder of files
        unlink(download_dir, recursive = TRUE, force = TRUE)
        message("DONE! Downloaded ", i, " genomes to ", combined_fasta)
    }
    return(combined_fasta)
}
