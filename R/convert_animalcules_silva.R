# Helper function
comp_line <- function(this_line, tax_cols) {
  new_line <- this_line
  empty_ind <- which(is.na(this_line))
  if (length(empty_ind) == length(tax_cols)) return(rep("unknown", length(tax_cols)))
  if (length(empty_ind) == 0) return(new_line)
  
  for (i in empty_ind) {
    filltax <- i - 1
    if (filltax %in% empty_ind) {
      new_line[i] <- new_line[filltax]
    } else {
      new_line[i] <- paste(tolower(substr(tax_cols[filltax], 1, 1)),
                           new_line[filltax], sep = "_")
    }
  }
  return(new_line)
}

# Helper function
read_in_id_silva <- function(path_id_counts, end_string,
                             species_headers_all) {
  name_file <- utils::tail(stringr::str_split(path_id_counts, "/")[[1]],
                           n = 1)
  output <- readr::read_csv(path_id_counts, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(.data$TaxonomyID)) %>%
    dplyr::select("read_count", "TaxonomyID", "Genome") %>%
    dplyr::mutate(sample = stringr::str_remove(name_file, end_string))
  counts_table <- output %>% dplyr::select(-"Genome")
  tax_table <- output %>% dplyr::select("TaxonomyID", "Genome")
    #species_headers_all %>%
    #dplyr::filter("TaxonomyID" %in% output$TaxonomyID)
  final_output <- list(counts_table, tax_table)
    return(final_output)
}

# Helper function
organize_tax_counts <- function(combined_list, tax_table) {
  order_ind <- match(combined_list$TaxonomyID, tax_table$TaxonomyID)
  tax_table_1 <- tax_table[order_ind, ]  %>%
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                function(x) stringr::str_remove_all(x, "\\[|\\]")),
                  dplyr::across(dplyr::everything(), stringr::str_to_sentence))
  tax_table_2 <- apply(tax_table_1, 1, comp_line,
                       tax_cols = colnames(tax_table_1)[seq_len(8)]) %>%
    t() %>%
    magrittr::set_rownames(tax_table_1$TaxonomyID) %>%
    magrittr::set_colnames(colnames(tax_table_1)) %>%
    as.data.frame()
  prelim_counts <- combined_list %>%
    tibble::column_to_rownames("TaxonomyID")
  # Upsample counts
  counts_up <- animalcules::upsample_counts(prelim_counts, tax_table_2, "species")
  tax_up <- tax_table_2[!duplicated(tax_table_2[, "species"]), ]
  rownames(tax_up) <- tax_up$species
  order_ind_2 <- match(rownames(counts_up), rownames(tax_up))
  tax_final <- tax_up[order_ind_2, ]
  return(list(counts = counts_up, tax = tax_final))
}

# Helper function
add_in_taxa_silva <- function(combined_pre, caching, path_to_write) {
  location <- "https://github.com/wejlab/metascope-docs/raw/main/all_silva_headers.rds"
  filename <- "all_silva_headers.rds"
  if (!caching) {
    if (!dir.exists(path_to_write)) dir.create(path_to_write)
    destination <- paste(path_to_write, filename, sep = "/")
    utils::download.file(location, destination)
  } else if (caching) {
    bfc <- .get_cache()
    rid <- BiocFileCache::bfcquery(bfc, filename, "rname")$rid
    if (!length(rid)) {
      rid <- names(BiocFileCache::bfcadd(bfc, filename, location))
    }
    if (!isFALSE(BiocFileCache::bfcneedsupdate(bfc, rid))) {
      BiocFileCache::bfcdownload(bfc, rid)
      BiocFileCache::bfcrpath(bfc, rids = rid)
    } else {
      message("Caching is set to TRUE, ",
              "and it appears that this file is already downloaded ",
              "in the cache. It will not be downloaded again.")
    }
    destination <- BiocFileCache::bfcrpath(bfc, rids = rid)
  }
  all_silva_headers <- readRDS(destination)
  tax_table_pre_1 <- lapply(combined_pre, function(x) x[[2]]) %>%
    dplyr::bind_rows()
  tax_table_pre <- tax_table_pre_1[!duplicated(tax_table_pre_1[, "TaxonomyID"]), ] %>%
    dplyr::left_join(all_silva_headers, by = c("TaxonomyID")) %>%
    dplyr::relocate("genus", "species", "TaxonomyID", .after = "family") %>%
    dplyr::select(-"Genome")
  return(tax_table_pre)
}

# Helper function
transform_cols <- function(x) {
  stringr::str_split_i(x, "__", i = -1) %>% stringr::str_replace_all("-", "_")
}

#' Create a multi-assay experiment from MetaScope output for usage with
#' animalcules with the SILVA 13_8 database
#'
#' Upon completion of the MetaScope pipeline, users can analyze and visualize
#' abundances in their samples using the animalcules package. This function
#' allows interoperability of \code{metascope_id} output with both animalcules
#' and QIIME. After running this function, the user should save the returned MAE
#' to an RDS file using a function like \code{saveRDS} to upload the output into
#' the \code{animalcules} package. NOTE: This function is for outputs that were
#' generated with the SILVA 13_8 database.
#'
#' @param meta_counts A vector of filepaths to the counts ID CSVs output by
#'   \code{metascope_id()} created with the SILVA database.
#' @inheritParams convert_animalcules
#' @inheritParams download_refseq
#' @returns Returns a MultiAssay Experiment file of combined sample counts data
#'   and/or saved biom file and mapping file for analysis with QIIME. The 
#'   MultiAssayExperiment will have a counts assay ("MGX").
#' @export
#' @importFrom rlang .data
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
#'   final_taxids <- c("AY846380.1.2583", "AY909584.1.2313", "HG531388.1.1375")
#'   final_genomes <- rep("Genome name", 3)
#'   best_hit <- sample(seq(100, 1050), 3)
#'   proportion <- best_hit/sum(best_hit) |> round(2)
#'   EMreads <- best_hit + round(runif(3), 1)
#'   EMprop <- proportion + 0.003
#'   dplyr::tibble("TaxonomyID" = final_taxids,
#'                 "Genome" = final_genomes,
#'                 "read_count" = best_hit, "Proportion" = proportion,
#'                 "EMreads" = EMreads, "EMProportion" = EMprop) |>
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
#' # Convert samples to MAE
#' outMAE <- convert_animalcules_silva(meta_counts = out_files,
#'                                     annot_path = annot_dat,
#'                                     which_annot_col = "Sample",
#'                                     end_string = ".metascope_id.csv",
#'                                     qiime_biom_out = FALSE,
#'                                     caching = TRUE)
#'
#' unlink(tempfolder, recursive = TRUE)
#' 

convert_animalcules_silva <- function(meta_counts, annot_path, which_annot_col,
                                      end_string = ".metascope_id.csv",
                                      qiime_biom_out = FALSE, path_to_write = ".",
                                      caching = TRUE) {
  combined_pre <- lapply(meta_counts, read_in_id_silva, end_string = end_string)
  combined_list <- lapply(combined_pre, function(x) x[[1]]) %>%
    data.table::rbindlist() %>% as.data.frame() %>%
    dplyr::mutate(sample = stringr::str_remove_all(sample, ".csv")) %>%
    dplyr::select("read_count", "TaxonomyID", "sample") %>%
    tidyr::pivot_wider(
      id_cols = "TaxonomyID", names_from = "sample",
      values_from = "read_count", values_fill = 0, id_expand = TRUE) %>%
    dplyr::group_by(.data$TaxonomyID) %>%
    dplyr::summarise(dplyr::across(dplyr::where(is.numeric), sum))
  
  tax_table_pre <- add_in_taxa_silva(combined_pre, caching, path_to_write)
  final_res <- organize_tax_counts(combined_list, tax_table_pre)
  MAE <- create_MAE(annot_path, which_annot_col, final_res$counts,
                    final_res$tax, path_to_write, qiime_biom_out)
  return(MAE)
}
