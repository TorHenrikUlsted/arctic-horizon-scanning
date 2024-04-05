filter_known = function(dfs, verbose = FALSE) {
  # The input parameters must be included and used
  # The return names must also be the same, make sure you return a list, even if it is only one return element
  # Everything else can be custom made, no need to us the help filters
  
  ##############
  # Initialize
  ##############
  
  # Change this to what you want to call your dataset. 
  # If you output the known species filter process as the same as known_name with an "_present" or "_absent" or both, then you do not need to initialize it yourself. 
  # Just make sure the filter process works like you want it to
  
  common_name <- "common"
  known_names <- c("known", "known2")
  
  dir_name <- paste0("./outputs/filter/", common_name)
  
  create_dir_if(dir_name)
  
  ###########################
  # Filter present & absent
  ###########################
  
  # Union_dfs merges and removes duplicates while also provide info on how many are removed
  
  combined_present <- data.table()
  combined_absent <- data.table()
  
  for (name in known_names) {
    present <- paste0(name, "_present")
    
    if (present %in% names(dfs)) {
      combined_present <- union_dfs(combined_present, dfs[[present]], verbose = verbose)
    }
    
    absent <- paste0(name, "_absent")
    
    if (absent %in% names(dfs))  {
      combined_absent <- union_dfs(combined_absent, dfs[[absent]], verbose = verbose)
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
    spec.out = paste0(data_common, "Absent"),
    fun = function() {
      
      combined_absent <-  dplyr::anti_join(combined_absent, combined_present, by = column)
      
      return(combined_absent)
    })
  
  return(list(
    present = combined_present,
    absent = combined_absent
  ))
}