filter_boreal <- function(spec.native, dfs, column, verbose = FALSE) {
  
  vebcat("Initiating boreal filter protocol.", color = "funInit")
  
  ##############
  # Initialize
  ##############
  
  arctic_present <- spec.native$present
  arctic_absent <- spec.native$absent
  glonaf_species <- dfs$glonaf_species
  boreal_species <- dfa$boreal_species
  
  boreal_dir <- "./outputs/filter/boreal"
  boreal_log <- paste0(boreal_dir, "/logs")
  
  create_dir_if(boreal_dir, boreal_log)
  
  # Here we remove the glonaf species, because we do not need duplicates
  boreal_species <- dplyr::anti_join(boreal_species, glonaf_species, by = column)
  
  
  ##################################################
  #           boreal arctic present absent         #
  ##################################################
  
  
  boreal_present <- write_filter_fun(
    file.out = paste0(boreal_dir, "/boreal-present-final.csv"),
    spec.in = boreal_species,
    spec.out = "Boreal present",
    fun = function() {
      # First merge to only get species from both dfs
      # This is because we do not want species not present in the arctic
      merged <- merge(boreal_species, arctic_present, by = column)
      
      return(merged)
    })
  
  # -------------------------------------------------------------------------- #
  
  boreal_absent <- write_filter_fun(
    file.out = paste0(boreal_dir, "/boreal-absent-final.csv"),
    spec.in = boreal_species,
    spec.out = "Boreal absent",
    fun = function() {
      # Remove species that are already present in the Arctic
      boreal_absent <-  dplyr::anti_join(boreal_species, arctic_present, by = column)
      
      # Update the absent list
      boreal_arc_absent = union_dfs(boreal_absent, arctic_absent, verbose = T)
      
      
      return(boreal_arc_absent)
    })
  
  vebcat("Boreal filter protocol completed successfully.", color = "funSuccess")
  
  return(list(
    spec = boreal_absent,
    dir = boreal_dir
  ))
}



