bt2_16S_params <- "-local -R 2 -N 0 -L 25 -i S,1,0.75 -k 10 --score-min L,100,1.28"
  # Code to prepare align_details list

# Save dataset
usethis::use_data(bt2_16S_params, overwrite = TRUE, internal = FALSE)
