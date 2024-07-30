filter_unknown <- function(known.filtered, dfs, column, verbose = FALSE) { # change the name of the function to 'filter_yournewname'

  ##############
  # Initialize
  ##############

  # Make sure the filter process works like you want it to. You can add custom functions inside the write_filter_fun.
  # Test the filter process by running it and checking the output files in the 'dir_name' directory.

  data_name <- "unknown" # Name the dataset whatever you want

  dir_name <- paste0("./outputs/filter/", data_name)

  create_dir_if(dir_name)

  mdwrite(
    config$files$post_seq_md,
    text = paste0("2;", data_name)
  )

  # Here you can get all species you have output from the filtering for the known species, change the names if necessary.

  known_present <- known.filtered$present
  known_absent <- known.filtered$absent

  data_species <- paste0(data_name, "_species")

  if (data_species %in% names(dfs)) {
    data_species <- dfs[[data_species]]
  } else {
    vebcat("Could not find the wrangled data.", color = "fatalError")
    catn("I am looking for", data_species, "in:")
    print(names(dfs))
  }


  ###########################
  # Filter present & absent
  ###########################

  # ------------------------------ Present ------------------------------- #
  # This is only used for information purposes. It is not used in the analysis. It generates a csv file.

  data_present <- write_filter_fun(
    file.out = paste0(dir_name, "/", data_name, "-present-final.csv"),
    spec.in = data_species,
    fun = function() {
      # First merge to only get species from both dfs
      data_present <- merge(data_species, known_present, by = column)

      return(data_present)
    }
  )


  # ------------------------------ Absent ------------------------------- #
  # This is the data used for the analysis.

  data_absent <- write_filter_fun(
    file.out = paste0(dir_name, "/", data_name, "-absent-final.csv"),
    spec.in = data_species,
    fun = function() {
      data_absent <- anti_union(data_species, known_present, by = column)

      # If you will use the known_absent in another list, remove it here instead
      data_absent <- union_dfs(data_absent, known_absent, by = column)

      return(data_absent)
    }
  )


  return(list(
    # These must remain 'spec' and 'dir'
    # 'data_absent' and 'dir_name' can be edited
    spec = data_absent,
    dir = dir_name
  ))
}
