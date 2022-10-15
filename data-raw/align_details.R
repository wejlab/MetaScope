# Code to prepare align_details list
align_details <- list(type = "dna", nthreads = 8, maxMismatches = 3,
                      nsubreads = 10, phredOffset = 33, unique = FALSE,
                      nBestLocations = 16)

# Save dataset
usethis::use_data(align_details, overwrite = TRUE, internal = FALSE)
