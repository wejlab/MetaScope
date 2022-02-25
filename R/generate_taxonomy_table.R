#' Helper function to create a taxonomy table
mk_table <- function (intable, taxon_ranks) {
  this_t <- as.data.frame(intable)
  if (!identical(this_t$rank, character(0))) {
    t_n <- nrow(this_t)
    if (this_t$rank[t_n] == "no rank" & this_t$rank[t_n - 1] == "species") {
      this_t$rank[t_n] <- "strain"
    }
    this_t %>%
      dplyr::filter(rank %in% taxon_ranks) %>% 
      dplyr::right_join(., dplyr::tibble(rank = taxon_ranks), 
                        by = "rank") %>% 
      dplyr::arrange(factor(rank, levels = taxon_ranks)) %>% 
      dplyr::select(name) %>% 
      .[seq_along(taxon_ranks), ]
  }
}

#' Generate the taxonomy table
#'
#' This code will generate the taxonomy table for all species or strains with
#' available NCBI genomes. Each column is a taxonomic rank or indication of a
#' strain; each row is a taxonomic classification for an unique species or
#' strain.
#' 
#' @importFrom magrittr %>%
#'
#' @return table of taxonomic relationships for all species or strains with
#' NCBI genomes.
#'
#' @examples
#' # Code not run
#' \dontrun{
#' generate_taxonomy_table()
#' }

generate_taxonomy_table <- function() {
  # Download the updated refseq table from NCBI
  # This table contains all species/strains with an available genome
  refseq_link <-
    "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
  tax_id <- read.table(refseq_link, header = T, sep = "\t",
                       comment.char = "", quote = "", skip = 1) %>%
    dplyr::distinct(taxid) %>% unlist %>% unname
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "strain")
  all_ncbi <- taxize::classification(tax_id, db = 'ncbi', max_tries = 3)
  taxonomy_table <- as.data.frame(t(dplyr::bind_rows(lapply(all_ncbi,
                                          mk_table,
                                          taxon_ranks))))
  colnames(taxonomy_table) <- taxon_ranks
  save(taxonomy_table, file = "data/taxonomy_table.rda")
}