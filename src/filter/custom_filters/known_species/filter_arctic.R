filter_arctic <- function(dfs, column, verbose = FALSE) {
  
  ##############
  # Initialize
  ##############
  
  mdwrite(
    post_seq_nums,
    heading = "2;Arctic"
  )
  
  aba_absent = dfs$aba_absent
  ambio_absent  = dfs$ambio_absent
  
  aba_present = dfs$aba_present
  ambio_present  = dfs$ambio_present
  
  # Union_dfs merges and removes duplicates while also provide info on how many are removed
  arctic_present <- union_dfs(aba_present, ambio_present, verbose = T)
  
  write_filter_fun(
    file.out = "./outputs/filter/arctic/arctic-present-final.csv",
    spec.in = arctic_present
  )
  
  # -------------------------------------------------------------------------- #
  
  arctic_absent <- union_dfs(aba_absent, ambio_absent, verbose = T)
  
  arctic_absent <- write_filter_fun(
    file.out = "./outputs/filter/arctic/arctic-absent-final.csv",
    spec.in = arctic_absent,
    fun = function() {
      
      # Also remove all arctic_present from absent in case some standard names have changed
      ab <- anti_union(arctic_absent, arctic_present, column)
      
      return(ab)
    })
  
  return(list(
    present = arctic_present,
    absent = arctic_absent
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
  
  write_filter_fun(
    file.out = paste0(dir_name, "/", common_name, "-present-final.csv"),
    spec.in = combined_present
  )
  
  # ---------------------------- Absent ----------------------------- #
  
  combined_absent <- write_filter_fun(
    file.out = paste0(dir_name, "/", common_name, "-absent-final.csv"),
    spec.in = combined_absent,
    fun = function() {
      
      combined_absent <-  anti_union(combined_absent, combined_present, column)
      
      return(combined_absent)
    })
  
  return(list(
    present = combined_present,
    absent = combined_absent
  ))
}
