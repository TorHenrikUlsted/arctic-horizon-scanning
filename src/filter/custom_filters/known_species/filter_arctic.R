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