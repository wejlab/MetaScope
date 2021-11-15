
#' Get child nodes from NCBI taxonomy
#'
#' This function will utilize a organism classification table to obtain all 
#' chidren species and/or strains with available NCBI reference sequences given
#' a parent taxon and its rank.
#'
#' @param input_taxon The parent taxon.
#' @param input_rank The taxonomic rank of the input taxon.
#' @param data A dataframe of orgnism classification information.
#' At minimum, should have a column indicating "strain", and and all others
#' should be taxonomic ranks. Each row should be a taxonomic relationship.
#' This defaults to the `taxonomy_table` object.
#'
#' @return Returns a vector of all the child species and/or strains of the
#' input taxon.
#'
#' @export
#'
#' @examples
#' ## Get all child species and strains in bacteria superkingdom
#' get_children('Bacteria','superkingdom')
#'
#' ## Get all child species and strains in fungi kingdom 
#' get_children('Fungi', 'kingdom')
#'
#' ## Get all child species in primate order 
#' get_children('Primates', 'order')
#'

get_children <- function(input_taxon, input_rank, data = taxonomy_table){
    # Get child strains
    ind <- tolower(data[, input_rank]) %in% tolower(input_taxon)
    strain_list <- unique(data[, "strain"][ind])
    strain_list <- strain_list[!is.na(strain_list)]
    # Get the children species
    ## Delete rows with child strains
    new_table <- data[!(data[, "strain"] %in% strain_list), ]
    ind <- tolower(new_table[, input_rank]) %in% tolower(input_taxon)
    species_list <- unique(new_table[, "species"][ind])
    species_list <- species_list[!is.na(species_list)]
    children_list <- c(strain_list,species_list)
    return(children_list)
}
