filter_arctic <- function(dfs, verbose = FALSE) {
  
  ##############
  # Initialize
  ##############
  
  aba_absent = dfs$aba_absent
  ambio_absent  = dfs$ambio_absent
  
  aba_present = dfs$aba_present
  ambio_present  = dfs$ambio_present
  
  # Union_dfs merges and removes duplicates while also provide info on how many are removed
  arctic_present <- union_dfs(aba_present, ambio_present, verbose = T)
  
  out_file <- "./outputs/filter/arctic/arctic-present-final.csv"
  create_dir_if(dirname(out_file))
  
  catn("Writing arctic_present to:", colcat(out_file, color = "output"))
  
  fwrite(arctic_present, out_file, bom = T)
  
  # -------------------------------------------------------------------------- #
  
  arctic_absent <- union_dfs(aba_absent, ambio_absent, verbose = T)
  
  arctic_absent <- write_filter_fun(
    file.out = "./outputs/filter/arctic/arctic-absent-final.csv",
    spec.in = arctic_absent,
    spec.out = "Arctic Absent",
    fun = function() {
      
      # Also remove all arctic_present from absent in case some standard names have changed
      ab <- dplyr::anti_join(arctic_absent, arctic_present, by = column)
      
      return(ab)
    })
  
  return(list(
    arctic_present = arctic_present,
    arctic_absent = arctic_absent
  ))
}

filter_arctictest = function(dfs, verbose = FALSE) {
  
  common_name <- "arctic"
  known_names <- c("aba", "ambio")
  
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