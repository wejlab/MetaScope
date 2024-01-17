# Create species headers data for silva

file_loc <- "/restricted/projectnb/pathoscope/work/tyler/MethodsAnalyses/DADA2/species_headers.txt"
species_headers <- readr::read_delim(file_loc,
                                     delim = " ", col_names = FALSE, show_col_types = FALSE) %>%
  dplyr::mutate(TaxonomyID = stringr::str_remove(X1, "^>")) %>%
  dplyr::relocate(TaxonomyID) %>%
  dplyr::mutate(Species = paste(X2, X3, sep = "_"), Genus = X2) %>%
  dplyr::select(TaxonomyID, Genus, Species) 
usethis::use_data(species_headers, overwrite = TRUE, internal = FALSE)
