filter_known <- function(dfs, verbose = FALSE) { # change the name of the function to 'filter_yournewname'

  ##############
  # Initialize
  ##############

  # Make sure the filter process works like you want it to.
  # Test the filter process by running it and checking the output files in the 'dir_name' directory.

  common_name <- "common_dataset_name" # Name the dataset whatever you want
  # If you have multiple datasets of known species, use the one below:
  known_names <- c("known_dataset", "known_dataset2")
  # If you only have one known dataset do:
  # known_names <- "dataset"

  dir_name <- paste0("./outputs/filter/", common_name) # This is where you can find the output csv files

  create_dir_if(dir_name)

  ###########################
  # Filter present & absent
  ###########################

  # Union_dfs and anti_union merges and removes duplicates while also providing info on how many are removed

  combined_present <- data.table()
  combined_absent <- data.table()

  for (name in known_names) {
    present <- paste0(name, "_present")

    if (present %in% names(dfs)) {
      combined_present <- union_dts(combined_present, dfs[[present]], verbose = verbose)
    }

    absent <- paste0(name, "_absent")

    if (absent %in% names(dfs)) {
      combined_absent <- union_dts(combined_absent, dfs[[absent]], verbose = verbose)
    }
  }

  # ---------------------------- Present ----------------------------- #

  out_file <- paste0(dir_name, "/", common_name, "-present-final.csv")

  catn("Writing", common_name, "_present to:", colcat(out_file, color = "output"))

  fwrite(combined_present, out_file, bom = T)

  # ---------------------------- Absent ----------------------------- #

  combined_absent <- write_filter_fun(
    file.out = paste0(dir_name, "/", common_name, "-absent-final.csv"),
    spec.in = combined_absent,
    fun = function() {
      combined_absent <- anti_union(combined_absent, combined_present, by = column)

      return(combined_absent)
    }
  )

  return(list(
    # These must remain 'present' and 'absent'
    # 'combined_present' and 'combined_absent' can be edited
    present = combined_present,
    absent = combined_absent
  ))
}
