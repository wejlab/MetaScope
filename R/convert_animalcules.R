# Helper functions for convert_animalcules

# Create biom file and mapping file
create_qiime_biom <- function(se_colData, taxonomy_table, which_annot_col,
                              counts_table, path_to_write) {
  # Consolidate counts based on genus
  counts_table_g <- counts_table %>%
    dplyr::mutate(genus = taxonomy_table$genus) %>%
    dplyr::group_by(.data$genus) %>%
    dplyr::summarise(dplyr::across(.fns = sum)) %>% as.data.frame()
  rownames(counts_table_g) <- counts_table_g$genus
  counts_table_g <- counts_table_g %>% dplyr::select(-.data$genus)
  # Match Indices
  ind <- base::match(se_colData[, which_annot_col], colnames(counts_table_g))
  counts_table_match <- counts_table_g[, c(ind)]
  # Formatting
  to_rm <- "[^[:alnum:] _.\\-+%;/:,]"
  func <- function(x) stringr::str_remove_all(x, to_rm)
  alt_col <- se_colData %>% dplyr::as_tibble() %>%
    dplyr::mutate("#SampleID" = se_colData[, which_annot_col],
                  # Sample Column: alphanumeric characters and periods ONLY
                  "#SampleID" = stringr::str_replace_all(.data$`#SampleID`, c(
                    "-" = "\\.", "_" = "\\.", " " = "\\.")),
                  "#SampleID" = stringr::str_remove_all(.data$`#SampleID`,
                                                        "[^[:alnum:].]")) %>%
    # Metadata columns: Only alphanumeric and [_.-+% ;:,/] characters
    apply(.data, 2, func) %>% ifelse(.data == " ", NA, .data) %>%
    ifelse(.data == "<NA>", NA, .data)
  # Remove unecessary columns
  ind <- colnames(alt_col) == which_annot_col
  alt_col2 <- alt_col[, -ind] %>% dplyr::as_tibble() %>%
    dplyr::mutate(BarcodeSequence = "-", LinkerPrimerSequence = "-",
                  Description = "-") %>% dplyr::relocate(
                    .data$`#SampleID`, .data$BarcodeSequence, .data$LinkerPrimerSequence)
  # Write files
  out_map <- paste(path_to_write, "QIIME_metadata_map.tsv", sep = "/")
  message("Writing TSV mapping file to ", out_map)
  utils::write.table(alt_col2, out_map, sep = "\t",
                     row.names = FALSE, quote = FALSE)
  # Write counts table to biom
  colnames(counts_table_match) <- alt_col2$`#SampleID`
  biom_obj <- biomformat::make_biom(counts_table_match)
  out_biom <- paste(path_to_write, "QIIME_featuretable.biom", sep = "/")
  message("Writing QIIME counts biom file to ", out_biom)
  biomformat::write_biom(biom_obj, out_biom)
}

# Create MAE Object
create_MAE <- function(annot_path, which_annot_col,
                       counts_table, taxonomy_table, path_to_write,
                       qiime_biom_out) {
  annot_dat <- readr::read_csv(annot_path, show_col_types = FALSE)
  se_colData <- annot_dat %>% # Only keep present samples in annotation data
    dplyr::mutate("sampcol" = unlist(annot_dat[, which_annot_col])) %>%
    dplyr::filter(.data$sampcol %in% colnames(counts_table)) %>%
    dplyr::select(-"sampcol") %>% S4Vectors::DataFrame()
  rownames(se_colData) <- se_colData[, which_annot_col]
  se_mgx <- counts_table %>% base::data.matrix() %>%
    S4Vectors::SimpleList() %>% magrittr::set_names("MGX")
  # Reorder colData according to se_mgx
  ind <- match(rownames(se_colData), colnames(se_mgx$MGX))
  se_colData <- se_colData[order(ind), ]
  se_rowData <- taxonomy_table %>% base::data.frame() %>%
    dplyr::mutate_all(as.character) %>% S4Vectors::DataFrame()
  microbe_se <- SummarizedExperiment::SummarizedExperiment(
    assays = se_mgx, colData = se_colData, rowData = se_rowData)
  MAE <- MultiAssayExperiment::MultiAssayExperiment(
    experiments = S4Vectors::SimpleList(MicrobeGenetics = microbe_se),
    colData = se_colData
  )
  if (qiime_biom_out) create_qiime_biom(
    se_colData, taxonomy_table, which_annot_col, counts_table, path_to_write)
  return(MAE)
}

# Read in the MetaScope_id CSVs
read_in_id <- function(path_id_counts, end_string, which_annot_col) {
  name_file <- utils::tail(stringr::str_split(path_id_counts, "/")[[1]],
                           n = 1)
  readr::read_csv(path_id_counts, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(.data$TaxonomyID)) %>%
    dplyr::select(.data$read_count, .data$TaxonomyID) %>%
    dplyr::mutate(sample = stringr::str_remove(name_file, end_string)) %>%
    return()
}

# Get input taxon phylogeny
class_taxon <- function(taxon, accession_path) {
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "strain")
  na_table <- data.frame(name = "Unknown", rank = taxon_ranks, id = 0)
  if (is.na(taxon)) return(na_table)

  classification_table <- taxonomizr::getTaxonomy(taxon,accession_path) |>
    as.data.frame()

  classification_table <- classification_table |>
    tidyr::pivot_longer(cols = colnames(classification_table), names_to = "rank", values_to = "name") |>
    dplyr::relocate(rank, .after = "name")
  classification_table$id <- c(rep(NA, nrow(classification_table) - 1), taxon)
  return(classification_table)
}

# Function written to account for accessions that weren't identified
# in metascope_id step
fill_in_missing <- function(combined_pre, accession_path) {
  na_ind <- combined_pre$TaxonomyID %>%
    as.double() %>% is.na() %>% which() %>% suppressWarnings()
  if(length(na_ind) > 0) {
    acc_list <- combined_pre$TaxonomyID[na_ind]
    result <- acc_list %>%
      find_accessions(quiet = TRUE, accession_path = accession_path) %>%
      plyr::aaply(1, function(x) x[1]) %>% unname()
    combined_pre$TaxonomyID[na_ind] <- result
  }
  return(combined_pre)
}

#' Create a multi-assay experiment from MetaScope output for usage with
#' animalcules
#'
#' Upon completion of the MetaScope pipeline, users can analyze and visualize
#' abundances in their samples using the animalcules package. This function
#' allows interoperability of \code{metascope_id} output with both animalcules
#' and QIIME. After running this function, the user should save the returned MAE
#' to an RDS file using a function like \code{saveRDS} to upload the output into
#' the \code{animalcules} package.
#'
#' @param meta_counts A vector of filepaths to the counts ID CSVs output by
#'   \code{metascope_id()}.
#' @param annot_path The filepath to the CSV annotation file for the samples.
#' This CSV metadata/annotation file should contain at least two columns,
#' one with names of all samples WITHOUT the extension listed in \code{end_string},
#' e.g. for output file "sample_x76.metascope_id.csv", the column specified in
#' \code{which_annot_col} should contain the entry "sample_x76". Sample names
#' containing characters "_", "-", and "." are fine, however sample names
#' beginning with numbers should be renamed to have a prefix, e.g. "777897sample"
#' should be renamed to "X777897sample" for both the output file name and the
#' annotation name.
#' @param end_string The end string used at the end of the metascope_id files.
#'   Default is ".metascope_id.csv".
#' @param which_annot_col The name of the column of the annotation file
#'   containing the sample IDs. These should be the same as the
#'   \code{meta_counts} root filenames.
#' @param qiime_biom_out Would you also like a qiime-compatible biom file
#'   output? If yes, two files will be saved: one is a biom file of the counts
#'   table, and the other is a specifically formatted mapping file of metadata
#'   information. Default is \code{FALSE}.
#' @param path_to_write If \code{qiime_biom_out = TRUE}, where should output QIIME
#'   files be written? Should be a character string of the folder path. Default is
#'   '.', i.e. the current working directory.
#' @param accession_path (character) Path to taxonomizr accessions. See
#'   \code{taxonomizr::prepareDatabase()}.
#' @returns Returns a MultiAssay Experiment file of combined sample counts data
#'   and/or biom file and mapping file for analysis with QIIME. The MultiAssay
#'   Experiment will have a counts assay ("MGX").
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' tempfolder <- tempfile()
#' dir.create(tempfolder)
#'
#' # Create three different samples
#' samp_names <- c("X123", "X456", "X789")
#' all_files <- file.path(tempfolder,
#'                        paste0(samp_names, ".csv"))
#'
#' create_IDcsv <- function (out_file) {
#'   final_taxids <- c("273036", "418127", "11234")
#'   final_genomes <- c(
#'     "Staphylococcus aureus RF122, complete sequence",
#'     "Staphylococcus aureus subsp. aureus Mu3, complete sequence",
#'     "Measles virus, complete genome")
#'   best_hit <- sample(seq(100, 1050), 3)
#'   proportion <- best_hit/sum(best_hit) |> round(2)
#'   EMreads <- best_hit + round(runif(3), 1)
#'   EMprop <- proportion + 0.003
#'   dplyr::tibble(TaxonomyID = final_taxids,
#'                 Genome = final_genomes,
#'                 read_count = best_hit, Proportion = proportion,
#'                 EMreads = EMreads, EMProportion = EMprop) |>
#'     dplyr::arrange(dplyr::desc(.data$read_count)) |>
#'     utils::write.csv(file = out_file, row.names = FALSE)
#'   message("Done!")
#'   return(out_file)
#' }
#' out_files <- vapply(all_files, create_IDcsv, FUN.VALUE = character(1))
#'
#' # Create annotation data for samples
#' annot_dat <- file.path(tempfolder, "annot.csv")
#' dplyr::tibble(Sample = samp_names, RSV = c("pos", "neg", "pos"),
#'               month = c("March", "July", "Aug"),
#'               yrsold = c(0.5, 0.6, 0.2)) |>
#'   utils::write.csv(file = annot_dat,
#'                    row.names = FALSE)
#'
#' # Create temporary taxonomizr accession
#' tmp_accession <- system.file("extdata", "example_accessions.sql", package = "MetaScope")
#'
#' # Convert samples to MAE
#' outMAE <- convert_animalcules(meta_counts = out_files,
#'                               annot_path = annot_dat,
#'                               which_annot_col = "Sample",
#'                               end_string = ".metascope_id.csv",
#'                               qiime_biom_out = FALSE,
#'                               accession_path = tmp_accession)
#'
#' unlink(tempfolder, recursive = TRUE)

convert_animalcules <- function(meta_counts, annot_path, which_annot_col,
                                end_string = ".metascope_id.csv",
                                qiime_biom_out = FALSE, path_to_write = ".",
                                accession_path = NULL) {
  combined_pre <- lapply(meta_counts, read_in_id, end_string = end_string,
                          which_annot_col = which_annot_col) %>%
    data.table::rbindlist() %>% dplyr::ungroup() %>% as.data.frame() %>%
    dplyr::mutate(sample = stringr::str_remove_all(sample, ".csv")) %>%
    dplyr::select("read_count", "TaxonomyID", "sample") %>%
    tidyr::pivot_wider(
      id_cols = .data$TaxonomyID, names_from = .data$sample,
      values_from = .data$read_count, values_fill = 0, id_expand = TRUE)
  # Which entries are not numeric, and try running them through genbank2uid
  combined_list <- fill_in_missing(combined_pre, accession_path) %>%
    dplyr::mutate("TaxonomyID" = as.numeric(.data$TaxonomyID))
  # Create taxonomy, counts tables
  taxon_ranks <- c("superkingdom", "kingdom", "phylum", "class", "order",
                   "family", "genus", "species", "strain")
  message("Looking up taxon UIDs in NCBI database")
  all_ncbi <- plyr::llply(combined_list$TaxonomyID, .fun = class_taxon,
                          accession_path = accession_path)
  # fix unknowns
  na_ind <- which(is.na(all_ncbi))
  unk_tab <- data.frame(name = "unknown", rank = taxon_ranks, id = 0)
  if (length(na_ind) > 0) for (i in na_ind) all_ncbi[[i]] <- unk_tab
  # Create table
  taxonomy_table <- plyr::llply(all_ncbi, mk_table, taxon_ranks) %>%
    dplyr::bind_rows() %>% as.data.frame()
  colnames(taxonomy_table) <- taxon_ranks
  counts_table <- combined_list %>% dplyr::select(-.data$TaxonomyID) %>%
    as.data.frame()
  # Remove any brackets
  taxonomy_table$species <- gsub("\\[|\\]", "", taxonomy_table$species)
  # Remove duplicated species
  if (sum(duplicated(taxonomy_table$species)) > 0) {
    ind <- which(duplicated(taxonomy_table$species))
    dup_sp <- taxonomy_table$species[ind]
    for (this_sp in dup_sp) {
      all_ind <- which(taxonomy_table$species == this_sp)
      counts_table[all_ind[1], ] <- base::colSums(counts_table[all_ind, ])
    }
    counts_table <- counts_table %>%
      dplyr::filter(!duplicated(taxonomy_table$species))
    taxonomy_table <- taxonomy_table %>%
      dplyr::filter(!duplicated(.data$species))
  }
  na_ind <- which(is.na(taxonomy_table$species))
  taxonomy_table$species[na_ind] <- paste0("g_", taxonomy_table$genus[na_ind])
  rownames(taxonomy_table) <- stringr::str_replace_all(taxonomy_table$species, " ", "_")
  rownames(counts_table) <- rownames(taxonomy_table)
  MAE <- create_MAE(
    annot_path, which_annot_col, counts_table,
    taxonomy_table, path_to_write, qiime_biom_out)
  return(MAE)
}
