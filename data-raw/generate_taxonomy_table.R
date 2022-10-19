find_taxonomy <- function(tax_id) {
  # Convert accessions to taxids and get genome names
  message("Obtaining taxonomy")
  # If URI length is greater than 2500 characters then split accession list
  URI_length <- nchar(paste(tax_id, collapse = "+"))
  if (URI_length > 2500) {
    chunks <- split(tax_id, ceiling(seq_along(tax_id) / 100))
    genome_return <- vector(mode = "list", length = length(chunks))
    message("Taxid list broken into ", length(chunks), " chunks")
    for (i in seq_along(chunks)) {
      success <- FALSE
      attempt <- 0
      # Attempt to get tax_id up to three times for each chunk
      while (!success) {
        try({
          attempt <- attempt + 1
          if (attempt > 1) message("Attempt #", attempt, " Chunk #", i)
          suppressMessages(
            tax_id_chunk <- taxize::classification(chunks[[i]], db = 'ncbi',
                                                   max_tries = 3))
          Sys.sleep(1)
          genome_return[[i]] <- tax_id_chunk
          success <- TRUE
        })
      }
    }
  } else suppressMessages(genome_return <- taxize::classification(
    taxid, db = 'ncbi', max_tries = 3))
  return(genome_return)
}

#' Generate the taxonomy table
#'
#' This code will generate the taxonomy table for all species or strains with
#' available NCBI genomes. Each column is a taxonomic rank or indication of a
#' strain; each row is a taxonomic classification for an unique species or
#' strain.
#' 
#' This returns a table of taxonomic relationships for all species or strains with
#' NCBI genomes.
#'
#' To use an NCBI key, add `ENTREZ_KEY = <your key here>` to your global
#' environment using `Sys.setenv()`

generate_taxonomy_table <- function() {
  # Download the updated refseq table from NCBI
  # This table contains all species/strains with an available genome
  refseq_link <-
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
  tax_id <- utils::read.table(refseq_link, header = TRUE, sep = "\t",
                              comment.char = "", quote = "", skip = 1) %>%
    dplyr::distinct(.data$taxid) %>% unlist() %>% unname()
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "strain")
  all_ncbi <- find_taxonomy(tax_id)
  taxonomy_table <- all_ncbi %>% unlist(recursive = FALSE) %>%
    lapply(mk_table, taxon_ranks = taxon_ranks) %>%
    dplyr::bind_rows() %>% t() %>% as.data.frame() %>%
    magrittr::set_colnames(taxon_ranks)
  usethis::use_data(taxonomy_table, internal = TRUE, overwrite = TRUE,
                    compress = "xz")
}

generate_taxonomy_table()
