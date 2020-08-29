#' Get the children taxons
#'
#' This function will utilize a organism classification table and get all 
#' the chidren species and/or strains with available NCBI refseq given 
#' a parent taxon and its rank.
#'
#' @param input_taxon The parent taxon.
#' @param input_rank The taxonomic rank of your input taxon.
#' @param data The dataframe of orgnism classification. 
#' Each column should be a taxnomic rank and each row should be a taxonomic relationship. 
#' Defaults to taxon_all table.
#' 
#' @return Returns a vector of all the children species and/or strains of the
#' input taxon 
#'
#' @examples
#' ## Get all the children species and strains under bacteria superkingdom.
#' get_children('Bacteria','superkingdom')
#'
#' ## Get all the children species and strains under fungi kingdom. 
#' get_children( 'Fungi', 'kingdom' )
#'
#' ## Get all the children species of primates order. 
#' get_children( 'Primates', 'order' )
#'
#' @export
#'

get_children <- function(input_taxon,
                         input_rank,
                         data){
  ## get the children strains
  strain_list <- unique(data[,"strain"][data[,input_rank] %in% input_taxon])
  strain_list <- strain_list[!is.na(strain_list)]
  ## get the children species
  ## delete the rows with children strains
  new_table <- data[!(data[,"strain"] %in% strain_list),]
  species_list <- unique(new_table[,"species"][new_table[,input_rank] %in% input_taxon])
  species_list <- species_list[!is.na(species_list)]
  children_list <- c(strain_list,species_list)
  return(children_list)
}
