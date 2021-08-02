
#' Generate the taxonomy table
#'
#' This code will generate the taxnonmy table for all species or strains with
#' available NCBI genomes. Each column is a taxonomic rank or indication of a
#' strain; each row is a taxonomic classification for an unique species or
#' strain.
#'
#' @return table of taxonomic relationships for all species or strains with
#' NCBI genomes.
#'
#' @examples
#' 
#' generate_taxonomy_table()
#' 


generate_taxonomy_table <- function() {
    # Download the updated refseq table from NCBI
    ## This table contains all species/strains with available genome
    refseq_link <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
    refseq_table <- utils::read.table(refseq_link, header = TRUE, sep = "\t", comment.char = "", quote = "", skip = 1)
    tax_id <- unique(refseq_table$taxid)
    n <- length(tax_id)
    taxonomy_table <- data.frame(matrix(NA, ncol = 9, nrow = length(tax_id)))
    taxon_ranks <- colnames(taxonomy_table) <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "strain")
    for (i in seq_len(n)) {
        this_id <- tax_id[i]
        if (i%%100 == 0) {
            print(paste("Obtained classification for ", i, " of ", n, ' species',
                        sep = ""))
        }
        suppressMessages(t <- taxize::classification(this_id, db = 'ncbi')[[1]])
        t_numrow <- nrow(t) 
        
        if (!identical(t$rank,character(0))) {
            ## Convert 'no rank' to 'strain'
            if (t$rank[t_numrow] == 'no rank' & t$rank[t_numrow - 1] == 'species') {
                t$rank[t_numrow] <- 'strain'
            }
            t <- t[which(t$rank %in% taxon_ranks), ]
            
            ## Find corresponding column positions for all parent taxons for given
            ## species/strain
            taxonomy_table[i, ] <- sapply(taxon_ranks, function(rn) ifelse(rn %in% t$rank, yes = subset(t$name, t$rank == rn), no = NA))
        }
    }
    taxonomy_table<- taxonomy_table[!apply(taxonomy_table, 1, function(x) all(is.na(x))), ]
    new <- taxonomy_table
    save(taxonomy_table, file = "data/taxonomy_table.rda")
}
