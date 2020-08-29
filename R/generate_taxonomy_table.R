#' Generate the taxonomy table
#' 
#' This code will generate the taxnonmy table for all species or strains with
#' available ncbi genomes. Each column is a taxonomic rank; each row is a taxonomic
#' classification for an unique species or strain.
#' 
#' 
require(taxize)
# download the updated refseq table from ncbi
refseq_link <- "ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt"
refseq_table <- read.table(refseq_link, header = T, sep = "\t",
                           comment.char = "",
                           quote = "", skip = 1)
# refseq_table contains all species/strains with available genome
# create an empty output table
taxonomy.table <- data.frame(superkingdom=character(),
                            kingdom=character(),
                            phylum=character(),
                            class=character(),
                            order=character(),
                            family=character(),
                            genus=character(),
                            species=character(),
                            strain=character(),
                            stringsAsFactors = F)

# filter all the unique taxid from refseq_table
tax_id <- unique(refseq_table$taxid)
for(i in tax_id){
  if(counter%%100 == 0){
    print(paste("Got classification for ",counter," out of ",length(tax_id), ' species',sep=""))
  }
  t <- classification(i,db='ncbi')[[1]]
  if(!identical(t$rank,character(0))){
    # convert 'no rank' to 'strain'
    if(t$rank[nrow(t)]=='no rank' & t$rank[nrow(t)-1]=='species'){
      t$rank[nrow(t)] <- 'strain'
    }
    t <- t[which(t$rank %in% taxon_ranks),]
    # a new table with the same columns
    new_table <- data.frame(superkingdom=character(),
                            kingdom=character(),
                            phylum=character(),
                            class=character(),
                            order=character(),
                            family=character(),
                            genus=character(),
                            species=character(),
                            strain=character(),
                            stringsAsFactors = F)
    new_row <- vector()
    # find the corresponding column positions for all the parent taxons for the species/strain
    for (j in c(1:length(colnames(new_table)))){
      rankname <- colnames(new_table)[j]
      if (rankname %in% t$rank){
        new_row <- c(new_row,t$name[t$rank==rankname])
      } else {new_row <- c(new_row,NA)}
    }
    # add this row to the new table
    new_table[1,] <- new_row
    #append the new table to the taxonomy table
    taxonomy.table <- rbind(taxonomy.table,new_table)
  }
  counter <- counter + 1
}