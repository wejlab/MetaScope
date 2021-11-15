
# Helper functions

#create_qiime_biom <- function() {
#  # Create mapping file
#  all_names <- c("#SampleID", "BarcodeSequence",	"LinkerPrimerSequence",
#                 colnames(sam_table), "Description")
#  biom_map <- data.frame(matrix("-", nrow = nrow(sam_table), ncol = length(all_names),
#                                dimnames = list(rownames(sam_table), all_names)))
#  biom_map[, 4:14] <- sam_table
#  biom_map$X.SampleID <- sam_table$Sample
#  biom_map <- biom_map[, -c(13)]
#  # write.csv(biom_map, "biom_map.csv")
#  colnames(biom_map)[1] <- "#SampleID"
#  
#  #biom_map %<>% 
#  #  mutate(shannon_div, invsimp_div, ginisim_div, unit_div)
#  
#  # biom_obj <- make_biom(counts_table, biom_map, tax_table)
#  biom_obj <- make_biom(counts_table)
#  write_biom(biom_obj, "QIIME_Work/exp2.biom")
#  write.table(biom_map, "QIIME_Work/biom_map.tsv", sep = "\t", row.names = FALSE)
#}

read_in_id <- function(path_id_counts, end_string, which_annot_col) {
  name_file <- utils::tail(stringr::str_split(path_id_counts, "/")[[1]], n = 1)
  meta_counts <- readr::read_csv(path_id_counts, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(TaxonomyID)) %>%
    dplyr::select(read_count, TaxonomyID) %>%
    dplyr::mutate(sample = stringr::str_remove(name_file, end_string))
}

#' Create a multi-assay experiment from MetaScope output for usage with animalcules

#' Upon completion of the MetaScope pipeline, users can analyze and visualize
#' abundances in their samples using the animalcules package. This function
#' allows interoperability of metascope_id output with both animalcules and
#' QIIME (coming soon!).

#' @param meta_counts A vector of filepaths to the counts ID CSVs output
#' by MetaScope
#' @param annot_path The filepath to the annotation file for the samples
#' @param end_string The end string used at the end of the metascope_id
#' files. Default is ".filtered.metascope_id.csv"
#' @param which_annot_col The name of the column of the annotation file
#' containing the sample IDs.
#' These should be the same as the \code{meta_counts} root filenames.
#' @param qiime_biom_out Would you also like a qiime-compatible biom file
#' output? Default is \code{FALSE}.
#' @returns returns a multi-assay experiment file of combined sample counts
#' data and/or biom file for QIIME. The multi-assay experiment will have
#' assays for the counts ("MGX"), log counts, CPM, and log CPM.
#' @export
#' @examples 
#' \donotrun{
#' all_files <- list.files("~/decamp/analysis/aodom/Novartis_COPD/MetaScope_run/Output", pattern = "*.filtered.metascope_id.csv", full.names = TRUE)
#' out <- convert_animalcules(meta_counts = all_files, annot_path = "~/decamp/analysis/aodom/Novartis_COPD/Data/DECAMP_patho_annotation.csv", end_string = ".filtered.metascope_id.csv", which_annot_col = "DECAMP_ID")
#' }
#'  # coming soon - this is from Aubrey
#'  

convert_animalcules <- function(meta_counts, annot_path, which_annot_col,
                                end_string = ".filtered.metascope_id.csv", 
                                qiime_biom_out = FALSE) {
  combined_list <- data.table::rbindlist(
    lapply(all_files, read_in_id, end_string = end_string,
           which_annot_col = which_annot_col)) %>%
    tidyr::pivot_wider(., id_cols = c("sample", "TaxonomyID"),
                       names_from = "sample", values_from = "read_count",
                       values_fill = 0)
  # Create taxonomy, counts tables
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "strain")
  all_ncbi <- taxize::classification(combined_list$TaxonomyID, db = "ncbi",
                                     max_tries = 4)
  taxonomy_table <- as.data.frame(t(sapply(all_ncbi, mk_table, taxon_ranks)))
  colnames(taxonomy_table) <- taxon_ranks
  counts_table <- combined_list %>% dplyr::select(-TaxonomyID) %>% 
    as.data.frame()
  # Remove duplicated species
  dup_sp <- taxonomy_table$species[which(duplicated(taxonomy_table$species))]
  all_ind <- which(taxonomy_table$species == dup_sp)
  counts_table[all_ind[1], ] <- base::colSums(counts_table[all_ind, ])
  counts_table %<>% dplyr::filter(!duplicated(taxonomy_table$species))
  taxonomy_table %<>% dplyr::filter(!duplicated(species))
  
  rownames(taxonomy_table) <- stringr::str_replace(taxonomy_table$species,
                                                   " ", "_")
  rownames(counts_table) <- rownames(taxonomy_table) 
  # Create MAE Object
  annot_dat <- readr::read_csv(annot_path, show_col_types = FALSE) 
  se_colData <- annot_dat %>% # Only keep present samples in annotation data
    dplyr::mutate(sampcol = unlist(annot_dat[, which_annot_col])) %>%
    dplyr::filter(sampcol %in% colnames(combined_list)) %>%
    dplyr::select(-sampcol) %>% S4Vectors::DataFrame()
  rownames(se_colData) <- se_colData[, which_annot_col]
  se_mgx <- counts_table %>% base::data.matrix() %>%
    S4Vectors::SimpleList() %>% magrittr::set_names("MGX")
  se_rowData <- taxonomy_table %>% base::data.frame() %>%
    dplyr::mutate_all(as.character) %>% S4Vectors::DataFrame()
  microbe_se <- SummarizedExperiment::SummarizedExperiment(
    assays = se_mgx, colData = se_colData, rowData = se_rowData) %>%
    TBSignatureProfiler::mkAssay(., input_name = "MGX", log = TRUE,
                                 output_name = "rawcounts")
  MAE <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(MicrobeGenetics = microbe_se),
    colData = se_colData)
}
