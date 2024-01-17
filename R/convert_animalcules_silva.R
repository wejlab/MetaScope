# Helper function
comp_line <- function(this_line, tax_cols) {
  empty_ind <- which(is.na(this_line))
  if (length(empty_ind) == length(tax_cols)) return(rep("unknown", length(tax_cols)))
  if (length(empty_ind) < 1) return(this_line)
  filltax <- min(empty_ind) - 1
  this_line[empty_ind] <- paste(tolower(substr(tax_cols[filltax], 1, 1)),
                                this_line[filltax], sep = "_")
  return(this_line)
}


# Helper function
read_in_id_silva <- function(path_id_counts, end_string) {
  name_file <- utils::tail(stringr::str_split(path_id_counts, "/")[[1]],
                           n = 1)
  output <- readr::read_csv(path_id_counts, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(.data$TaxonomyID)) %>%
    dplyr::select(.data$read_count, .data$TaxonomyID, .data$Genome) %>%
    dplyr::mutate(sample = stringr::str_remove(name_file, end_string))
  split_taxa <- stringr::str_split(output$Genome, pattern = ";")
  domain_names <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
  # Grab first 6 elements of list
  up_to_genus <- lapply(split_taxa, function(x) x[1:length(domain_names)])
  counts_table <- output %>% dplyr::select(-Genome)
  tax_table <- up_to_genus %>%
    unlist() %>% matrix(ncol = length(domain_names), byrow = TRUE) %>%
    data.frame() %>%
    magrittr::set_colnames(domain_names) %>%
    dplyr::mutate(TaxonomyID = counts_table$TaxonomyID) %>%
    tibble::as_tibble()
  final_output <- list(counts_table, tax_table)
    return(final_output)
}

# Helper function
organize_tax_counts <- function(combined_list, tax_table) {
  
  order_ind <- match(combined_list$TaxonomyID, tax_table$TaxonomyID)
  tax_table_new <- tax_table[order_ind, ]
  
  tax_table_3 <- tax_table_new %>%
    dplyr::select(-TaxonomyID) %>%
    dplyr::mutate(Species = stringr::str_remove(Species, pattern = "^s__"),
                  Species = stringr::str_replace(Species, pattern = "_", "_"),
                  Species = stringr::str_replace(Species, pattern = "_", "_"),
                  Species = stringr::str_replace(Species, pattern = " ", "_"))
  
  tax_table_4 <- apply(tax_table_3, 2,
                       function(x) stringr::str_split_i(x, "__", i = -1)) %>%
    as.data.frame()
  tax_cols <- colnames(tax_table_4)
  tax_table_5 <- apply(tax_table_4, 1, comp_line, tax_cols = tax_cols) %>% t() %>%
    magrittr::set_rownames(tax_table_new$TaxonomyID) %>%
    magrittr::set_colnames(tax_cols) %>%
    as.data.frame()
  counts_final <- combined_list %>%
    dplyr::mutate(Species = tax_table_5$Species) %>%
    dplyr::select(-TaxonomyID) %>%
    dplyr::group_by(Species) %>%
    dplyr::summarise(across(where(is.numeric), sum)) %>%
    as.data.frame() %>%
    dplyr::arrange(as.character(Species)) #%>%
    #dplyr::rename_with(transform_cols)
  rownames(counts_final) <- counts_final$Species
  counts_final <- dplyr::select(counts_final, -Species)
  tax_final <- tax_table_5 %>%
    dplyr::distinct(Species, .keep_all = TRUE) %>%
    dplyr::arrange(Species)
  rownames(tax_final) <- tax_final$Species
  return(list(counts = counts_final, tax = tax_final))
}

# Helper function
transform_cols <- function(x) {
  stringr::str_split_i(x, "__", i = -1) %>% stringr::str_replace_all("-", "_")
}

#' Create a multi-assay experiment from MetaScope output for usage with
#' animalcules with the SILVA database
#'
#' Upon completion of the MetaScope pipeline, users can analyze and visualize
#' abundances in their samples using the animalcules package. This function
#' allows interoperability of \code{metascope_id} output with both animalcules
#' and QIIME. After running this function, the user should save the returned MAE
#' to an RDS file using a function like \code{saveRDS} to upload the output into
#' the \code{animalcules} package. NOTE: This function is for outputs that were
#' generated with the SILVA database.
#'
#' @param meta_counts A vector of filepaths to the counts ID CSVs output by
#'   \code{metascope_id()} created with the SILVA database.
#' @inheritParams convert_animalcules
#' @returns Returns a MultiAssay Experiment file of combined sample counts data
#'   and/or biom file and mapping file for analysis with QIIME. The MultiAssay
#'   Experiment will have a counts assay ("MGX").
#' @export
#' @importFrom rlang .data
#' 
convert_animalcules_silva <- function(meta_counts, annot_path, which_annot_col,
                                      end_string = ".metascope_id.csv",
                                      qiime_biom_out = FALSE, path_to_write = ".") {
  combined_pre <- lapply(meta_counts, read_in_id_silva, end_string = end_string)
  combined_list <- lapply(combined_pre, function(x) x[[1]]) %>%
    data.table::rbindlist() %>% as.data.frame() %>%
    dplyr::mutate(sample = stringr::str_remove_all(sample, ".csv")) %>%
    dplyr::select("read_count", "TaxonomyID", "sample") %>%
    tidyr::pivot_wider(
      id_cols = .data$TaxonomyID, names_from = .data$sample,
      values_from = .data$read_count, values_fill = 0, id_expand = TRUE) %>%
    dplyr::group_by(TaxonomyID) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum))
  data_env <- new.env(parent = emptyenv())
  utils::data("species_headers", envir = data_env, package = "MetaScope")
  species_headers_all <- data_env[["species_headers"]] %>%
    dplyr::filter(TaxonomyID %in% combined_list$TaxonomyID)
  tax_table_pre <- lapply(combined_pre, function(x) x[[2]]) %>%
    data.table::rbindlist() %>%
    dplyr::distinct(TaxonomyID, .keep_all = TRUE) %>%
    dplyr::left_join(species_headers_all, by = "TaxonomyID") %>%
    dplyr::mutate(Genus = ifelse(!is.na(Genus.y), Genus.y, Genus.x)) %>%
    dplyr::select(-c(Genus.y, Genus.x)) %>%
    dplyr::relocate(Species, .after = Genus)
  
  final_res <- organize_tax_counts(combined_list, tax_table_pre)
  counts_table <- final_res$counts
  taxonomy_table <- final_res$tax
  
  MAE <- create_MAE(annot_path, which_annot_col, counts_table,
    taxonomy_table, path_to_write, qiime_biom_out)
  return(MAE)
}
