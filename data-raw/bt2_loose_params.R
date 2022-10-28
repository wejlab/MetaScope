bt2_loose_params <- "--local -k 100 -D 20 -R 3 -L 3 -N 1 -p 8 --gbar 1 --mp 3"
# Code to prepare align_details list

# Save dataset
usethis::use_data(bt2_loose_params, overwrite = TRUE, internal = FALSE)
