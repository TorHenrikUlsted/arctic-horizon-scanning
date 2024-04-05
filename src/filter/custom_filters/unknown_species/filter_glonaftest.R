filter_glonaftest = function(known.filtered, dfs, verbose = FALSE) {
  ##############
  # Initialize
  ##############
  
  # Change this to what you want to call your dataset. 
  # Just make sure the filter process works like you want it to you can add custom functions inside the write_filter_fun
  
  data_name <- "glonaf"
  
  dir_name <- paste0("./outputs/filter/", data_name)
  
  create_dir_if(dir_name)
  
  # Here you can get all species you have output from the filtering for the known species, change the names if necessary.
  
  known_present <- known$present
  known_absent <- known$absent
  
  data_species <- paste0(data_name, "_species")
  
  if (data_species %in% names(dfs)) {
    data_species <- dfs[[data_species]]
  } 
  
  ###########################
  # Filter present & absent
  ###########################
  
  # ------------------------------ Present ------------------------------- #
  
  data_present <- write_filter_fun(
    file.out = paste0(dir_name, "/", data_name, "-present-final.csv"),
    spec.in = data_species,
    spec.out = paste0(data_name, "Present"),
    fun = function() {
      # First merge to only get species from both dfs
      data_present <- merge(data_species, known_present, by = column)
      
      return(data_present)
    })
  
  
  # ------------------------------ Absent ------------------------------- #
  
  data_absent <- write_filter_fun(
    file.out = paste0(dir_name, "/", data_name, "-absent-final.csv"),
    spec.in = data_species,
    spec.out = paste0(data_name, "Absent"),
    fun = function() {
      data_absent <-  dplyr::anti_join(data_species, known_present, by = column)
      
      # Remove species in the absent list as well, and move those to the boreal species
      data_absent <-  dplyr::anti_join(data_absent, known_absent, by = column)
      
      return(data_absent)
    })
  
  
  return(list(
    spec = data_absent,
    dir = dir_name
  ))
}