bt2_missing_params <- "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.4"
# Code to prepare align_details list

# Save dataset
usethis::use_data(bt2_missing_params, overwrite = TRUE, internal = FALSE)
