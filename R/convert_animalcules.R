
# Helper functions for convert_animalcules

# Create biom file and mapping file
  # TO ADD: allow genus/other level combining???
create_qiime_biom <- function(se_colData, taxonomy_table, which_annot_col,
                              counts_table, path_to_write) {
  # Create mapping file
  all_names <- c("#SampleID", "BarcodeSequence",	"LinkerPrimerSequence",
                 colnames(se_colData), "Description")
  biom_map <- data.frame(matrix("-", nrow = nrow(se_colData),
                                ncol = length(all_names),
                                dimnames = list(c(), c())))
  colnames(biom_map) <- all_names
  biom_map[, colnames(se_colData)] <- as.matrix(se_colData)
  # Formatting
  alt_col <- se_colData %>% dplyr::as_tibble() %>%
    dplyr::mutate(`#SampleID` = se_colData[, which_annot_col],
           # Sample Column: alphanumeric characters and periods ONLY
           `#SampleID` = stringr::str_replace_all(
      `#SampleID`, c("-" = "\\.", "_" = "\\.", " " = "\\.")),
      `#SampleID` = stringr::str_remove_all(`#SampleID`, "[^[:alnum:].]")) %>%
    # Metadata columns: Only alphanumeric and [_.-+% ;:,/] characters
    apply(., 2, function(x) stringr::str_remove_all(
      x, "[^[:alnum:] _.\\-+%;/:,]")) %>%
    ifelse(. == " ", NA, .) # Change " " to NA
  # Remove unecessary columns
  ind <- colnames(se_colData) == which_annot_col
  alt_col2 <- alt_col[, -ind] %>% dplyr::as_tibble() %>% 
    dplyr::mutate(BarcodeSequence = "-", LinkerPrimerSequence = "-",
                  Description = "-") %>%
    dplyr::relocate(`#SampleID`, BarcodeSequence, LinkerPrimerSequence)
  # Write files
  out_map <- paste(path_to_write, "QIIME_metadata_map.tsv", sep = "/")
  message("Writing mapping TSV to ", out_map)
  write.table(biom_map, out_map, sep = "\t", row.names = FALSE, quote = FALSE)
  biom_obj <- biomformat::make_biom(counts_table)
  out_biom <- paste(path_to_write, "QIIME_featuretable.biom", sep = "/")
  message("Writing QIIME counts biom file to ", out_biom)
  biomformat::write_biom(biom_obj, out_biom)
}

# Create MAE Object
create_MAE <- function(annot_path, which_annot_col, combined_list,
                       counts_table, taxonomy_table, path_to_write,
                       qiime_biom_out) {
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
  out_MAE <- paste(path_to_write, "metascope_animalcules_MAE.RDS", sep = "/")
  message("Writing animalcules Multi-Assay Experiment to ", out_MAE)
  saveRDS(MAE, out_MAE)
  if (qiime_biom_out) create_qiime_biom(se_colData, taxonomy_table,
                                        which_annot_col, counts_table,
                                        path_to_write)
  return(MAE)
}

# Read in the MetaScope_id CSVs
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
#' output? If yes, two files will be saved: one is a boDefault is \code{FALSE}.
#' @param path_to_write Where should output animalcules and/or QIIME files
#' be written to? Should be a character string of the folder path.
#' Default is '.', i.e. the current working directory.
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
                                qiime_biom_out = FALSE,
                                path_to_write = ".") {
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
  if (sum(duplicated(taxonomy_table$species)) > 0) {
    dup_sp <- taxonomy_table$species[which(duplicated(taxonomy_table$species))]
    for (this_sp in dup_sp) {
      all_ind <- which(taxonomy_table$species == this_sp)
      counts_table[all_ind[1], ] <- base::colSums(counts_table[all_ind, ])
    }
    counts_table %<>% dplyr::filter(!duplicated(taxonomy_table$species))
    taxonomy_table %<>% dplyr::filter(!duplicated(species))
  }
  rownames(taxonomy_table) <- stringr::str_replace(taxonomy_table$species,
                                                   " ", "_")
  rownames(counts_table) <- rownames(taxonomy_table) 
  MAE <- create_MAE(annot_path, which_annot_col, combined_list, counts_table,
             taxonomy_table, path_to_write, qiime_biom_out)
  return(MAE)
}
