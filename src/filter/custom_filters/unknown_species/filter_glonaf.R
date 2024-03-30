filter_glonaf = function(spec.known, dfs, verbose = FALSE) {
  ##############
  # Initialize
  ##############
  glonaf_species <- dfs$glonaf_species
  glonaf_dir <- "./outputs/filter/glonaf"
  
  create_dir_if(glonaf_dir)
  
  
  ##############
  # Filter
  ##############
  
  glonaf_present <- write_filter_fun(
    file.out = paste0(glonaf_dir, "/glonaf-present-final.csv"),
    spec.in = glonaf_species,
    spec.out = "Glonaf Present",
    fun = function() {
      # First merge to only get species from both dfs
      glonaf_present <- merge(glonaf_species, spec.known$present, by = column)
      
      return(glonaf_present)
    })
  
  
  # -------------------------------------------------------------------------- #
  
  glonaf_absent <- write_filter_fun(
    file.out = paste0(glonaf_dir, "/glonaf-absent-final.csv"),
    spec.in = glonaf_species,
    spec.out = "Glonaf Absent",
    fun = function() {
      glonaf_absent <-  dplyr::anti_join(glonaf_species, spec.known$present, by = column)
      
      # Remove species in the absent list as well, and move those to the boreal species
      glonaf_absent <-  dplyr::anti_join(glonaf_species, spec.known$absent, by = column)
      
      return(glonaf_absent)
    })
  
  return(list(
    spec = glonaf_absent,
    dir = glonaf_dir
  ))
}