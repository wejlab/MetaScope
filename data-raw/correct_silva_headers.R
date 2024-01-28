correct_silva_headers <- function(all_tax) {
  taxonomy_table <- get0("taxonomy_table", envir = asNamespace("MetaScope"))
  
  # Split all headers by ";" delimiter
  all_tax_split <- all_tax$Genome |>
    plyr::aaply(1, function(x) stringr::str_split(x, ";")) |>
    magrittr::set_names(all_tax$TaxonomyID) |>
    plyr::llply(function(x) gsub("Proteobacteria", "Pseudomonadota", x)) |>
    plyr::llply(function(x) gsub("Cyanobacteria", "Cyanobacteriota", x))
  
  # for each element, check where it belongs in taxonomy_table
  # Do this via a series of vectorized functions
  in_fun <- function(x) {
    out <- unname(unlist(unique(x)))
    return(out[!is.na(out)])
  }
  all_uniq_tax <- plyr::alply(taxonomy_table, 2, in_fun, .parallel = FALSE) |>
    magrittr::set_names(colnames(taxonomy_table))
  
  find_element <- function(element, search_string) {
    if (search_string %in% all_uniq_tax[[element]]) {
      return(names(all_uniq_tax)[element])
    } else {
      return(NA)
    }
  }
  
  get_tax_level <- function(silva_tax_level) {
    result <- plyr::ldply(seq_along(all_uniq_tax),
                          find_element, search_string = silva_tax_level,
                          .parallel = FALSE) |>
      unlist()
    if (!all(is.na(result))) return(result[!is.na(result)][1])
    return(NA)
  }
  
  grab_taxa_overhead <- function(genome_name) {
    out <- tibble::tibble(input = genome_name) |>
      plyr::adply(1, get_tax_level, .id = "input", .parallel = FALSE)
    return(out)
  }
  final_desig <- plyr::llply(all_tax_split, grab_taxa_overhead, .parallel = FALSE)
  all_taxa_cols <- c("superkingdom", "kingdom",
                     "phylum", "class", "order", "family", "genus")
  all_taxa_bact <- c("superkingdom", "phylum", "class", "order", "family", "genus")
  
  format_taxa <- function(slice) {
    slice <- slice |>
      dplyr::filter(stringr::str_detect(input, "[[:alnum:]]"))
    ind <- is.na(slice[, 2])
    
    if (slice[1, 1] == "Bacteria") {
      which_taxa <- all_taxa_bact
    } else which_taxa <- all_taxa_cols
    
    # If it's sandwiched between valid taxa, make it right
    for (k in seq_along(ind)) {
      if (isFALSE(ind[k])){
        NULL 
      } else if (k > 1) {
        check <- isFALSE(ind[k - 1]) && isFALSE(ind[k + 1])
        if (check && k != length(ind)) {
          left <- which(slice[k - 1, 2] == which_taxa)
          right <- which(slice[k + 1, 2] == which_taxa)
          if (right - left == 2) {
            this_in <- which_taxa[left + 1]
            if (this_in %in% slice[, 2]) {
              slice[k, 2] <- NA
            } else slice[k, 2] <- this_in
          }
        } else {
          if (isFALSE(ind[k - 1])) {
            is_us <- stringr::str_detect(slice[k, 1], "_")
            is_un <- stringr::str_detect(slice[k, 1], "uncultured")
            is_num <- stringr::str_detect(slice[k, 1], "\\d+")
            is_sp <-  stringr::str_detect(slice[k, 1], "sp.")
            left <- which(slice[k - 1, 2] == which_taxa)
            if (!is_us && !is_un && !is_num && !is_sp) {
              this_in <- which_taxa[left + 1]
              if (this_in %in% slice[, 2]) {
                slice[k, 2] <- NA
              } else slice[k, 2] <- this_in
            }
          }
        }
      }
      ind <- is.na(slice[, 2])
    }
    slice_out <- slice %>% dplyr::filter(!is.na(.data$V11)) %>%
      dplyr::distinct(.data$V11, .keep_all = TRUE)
    return(slice_out)
  }
  pre_out <- plyr::llply(final_desig, format_taxa) |>
    dplyr::bind_rows(.id = "taxid") |>
    tidyr::pivot_wider(id_cols = "taxid", names_from = "V11",
                       values_from = "input") |>
    dplyr::relocate(dplyr::all_of(all_taxa_cols), "taxid") |>
    dplyr::rename("TaxonomyID" = "taxid")
  
  # For parallel
  return(pre_out)
}

file_loc <- "/restricted/projectnb/pathoscope/reflib/2019_silva/all_headers.txt"
silva_headers <- readr::read_delim(file_loc,
                                     delim = "\t", col_names = FALSE, show_col_types = FALSE) |>
  dplyr::mutate(X1 = stringr::str_replace(X1, ";", "@")) |>
  tidyr::separate(X1, into = c("TaxonomyID", "Genome"), sep = "@")

# How long
510984 # total taxa
# 50 taxa / sec
# so ~ 170 minutes

now <- Sys.time()
silva_headers_final <- correct_silva_headers(silva_headers)
Sys.time() - now

# Add in silva species headers
file_loc <- "/restricted/projectnb/pathoscope/work/tyler/MethodsAnalyses/DADA2/species_headers.txt"
species_headers <- readr::read_delim(file_loc,
                                     delim = " ", col_names = FALSE, show_col_types = FALSE) |>
  dplyr::mutate(TaxonomyID = stringr::str_remove(X1, "^>")) |>
  dplyr::relocate(TaxonomyID) |>
  dplyr::mutate(species = paste(X2, X3, sep = "_"), genus = X2) |>
  dplyr::select(TaxonomyID, genus, species)

all_silva_headers <- dplyr::left_join(silva_headers_final, species_headers, by = c("TaxonomyID")) |>
  dplyr::mutate("genus" = replace(.data$genus.x, !is.na(.data$genus.y),
                                  .data$genus.y[!is.na(.data$genus.y)])) |>
  dplyr::select(!tidyselect::starts_with("genus.")) |>
  dplyr::relocate("genus", "species", "TaxonomyID", .after = "family")

# Too big so we save to docs
#usethis::use_data(all_silva_headers, overwrite = TRUE, internal = FALSE,
#                  compress = "bzip2")
saveRDS(all_silva_headers, file = "~/pathoscope/work/aubrey/MetaScope/Meta/all_silva_headers.rds",
        compress = "bzip2")