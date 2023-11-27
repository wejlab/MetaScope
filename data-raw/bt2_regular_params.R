bt2_regular_params <- "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.7"
# Code to prepare align_details list

# Save dataset
usethis::use_data(bt2_regular_params, overwrite = TRUE, internal = FALSE)
