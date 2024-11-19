#' Create a multi-assay experiment from PathoScope 2.0 output for usage with
#' animalcules
#'
#' This function serves as a legacy integration method for usage with
#' PathoScope 2.0 outputs. Upon completion of the PathoScope 2.0 pipeline, users can analyze and visualize
#' abundances in their samples using the animalcules package. After running this function, the user should save the returned MAE
#' to an RDS file using a function like \code{saveRDS} to upload the output into
#' the \code{animalcules} package.
#'
#' @param patho_counts Character string, a directory filepath to the counts ID CSVs output by
#'   \code{metascope_id()}.
#' @param annot_path The filepath to the CSV annotation file for the samples.
#' @param end_string The end string used at the end of the metascope_id files.
#'   Default is ".metascope_id.csv".
#' @param which_annot_col The name of the column of the annotation file
#'   containing the sample IDs. These should be the same as the
#'   \code{meta_counts} root filenames.
#' @returns Returns a MultiAssay Experiment file of combined sample counts data.
#' The MultiAssay Experiment will have a counts assay ("MGX").
#' @export
#' @importFrom rlang .data
#' 
convert_animalcules_patho <- function(patho_counts, annot_path,
                                      which_annot_col,
                                      end_string = "-sam-report.tsv") {
  count_table <- animalcules::read_pathoscope_data(
    input_dir = patho_counts,
    pathoreport_file_suffix = end_string)$countdat
  # Choose only the samples in metadata that have counts data as well
  metadata_table <- readr::read_csv(annot_path, show_col_types = FALSE) |>
    as.data.frame()
  rownames(metadata_table) <- metadata_table[, which_annot_col]
  sample_overlap <- base::intersect(colnames(count_table), rownames(metadata_table))
  if (length(sample_overlap) < length(colnames(count_table))) {
    no_info <- paste(colnames(count_table)[which(!colnames(count_table) %in% sample_overlap)],
                     collapse = ",")
    message("The following samples don't have metadata info:", no_info)
    
    count_table <- count_table[,which(colnames(count_table) %in% sample_overlap)]
  }
  metadata_table <- metadata_table[match(colnames(count_table),
                                         rownames(metadata_table)), , drop = FALSE]
  # Test and fix the constant/zero row
  row.remove.index <- c()
  if (sum(base::rowSums(as.matrix(count_table)) == 0) > 0){
    row.remove.index <- which(base::rowSums(as.matrix(count_table)) == 0)
    count_table <- count_table[-row.remove.index,]
  }
  ids <- rownames(count_table)
  tids <- unlist(lapply(ids, FUN = animalcules::grep_tid))
  if (sum(is.na(tids)) > 0){
    tid_remove <- which(is.na(tids))
    ids <- ids[-tid_remove]
    tids <- tids[-tid_remove]
    count_table <- count_table[-tid_remove,]
  }
  taxonLevels <- animalcules::find_taxonomy(tids)
  tax_table <- animalcules::find_taxon_mat(ids, taxonLevels)
  # replace spaces in tax name with underscore
  tax_table <- as.data.frame(apply(tax_table,
                                   2,
                                   function(x) gsub('\\s+', '_',x)))
  
  # create MAE object
  se_mgx <-
    count_table |>
    base::data.matrix() |>
    S4Vectors::SimpleList() |>
    magrittr::set_names("MGX")
  
  # Read in metadata
  se_colData <-
    metadata_table |>
    S4Vectors::DataFrame()
  
  se_rowData <-
    tax_table |>
    base::data.frame() |>
    dplyr::mutate_all(as.character) |>
    dplyr::select("superkingdom", "phylum", "class", "order",
                  "family", "genus", "species") |>
    S4Vectors::DataFrame()
  
  microbe_se <-
    SummarizedExperiment::SummarizedExperiment(assays = se_mgx,
                                               colData = se_colData,
                                               rowData = se_rowData)
  mae_experiments <-
    S4Vectors::SimpleList(MicrobeGenetics = microbe_se)
  
  MAE <-
    MultiAssayExperiment::MultiAssayExperiment(experiments = mae_experiments,
                                               colData = se_colData)
  
  return(MAE)
}
