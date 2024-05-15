filter_glonaf = function(known.filtered, dts, column, verbose = FALSE) {
  ##############
  # Initialize
  ##############
  mdwrite(
    post_seq_nums,
    heading = "2;GloNAF"
  )
  
  glonaf_species <- dts$glonaf_absent
  glonaf_dir <- "./outputs/filter/glonaf"
  create_dir_if(glonaf_dir)
  
  arctic_present <- known.filtered$present
  arctic_absent <- known.filtered$absent
  
  mdwrite(
    post_seq_nums,
    heading = "2;GloNAF"
  )
  
  ##############
  # Filter
  ##############
  
  glonaf_present <- write_filter_fun(
    file.out = paste0(glonaf_dir, "/glonaf-present-final.csv"),
    spec.in = glonaf_species,
    fun = function() {
      # First merge to only get species from both dts
      glonaf_present <- union_dfs(glonaf_species, arctic_present)
      
      return(glonaf_present)
    })
  
  
  # -------------------------------------------------------------------------- #
  
  glonaf_absent <- write_filter_fun(
    file.out = paste0(glonaf_dir, "/glonaf-absent-final.csv"),
    spec.in = glonaf_species,
    fun = function() {
      glonaf_absent <-  anti_union(glonaf_species, arctic_present, column)
      
      # Remove species in the absent list as well, and move those to the boreal species
      glonaf_absent <-  anti_union(glonaf_absent, arctic_absent, column)
      
      return(glonaf_absent)
    })
  
  return(list(
    spec = glonaf_absent,
    dir = glonaf_dir
  ))
}