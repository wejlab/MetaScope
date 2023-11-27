# Code to prepare align_details list
bt2_16S_params <-
  "--local -R 2 -N 0 -L 25 -i S,1,0.75 -k 5 --score-min L,0,1.88"

# Save dataset
usethis::use_data(bt2_16S_params, overwrite = TRUE, internal = FALSE)
