analyze_hv_stats <- function(region_hv, sp_hv, spec.name, verbose) {
  cat("Analyzing hypervolume for", cc$lightSteelBlue(sp_hv@Name), "species. \n")
  
  hv_set <- hypervolume_set(sp_hv, region_hv, check.memory = F)
  
  hv_stats <- hypervolume_overlap_statistics(hv_set)
  
  sp_surviv_region <- 1 - hv_stats[[4]]
  
  cat("Volume of", spec.name, "overlap in the CAVM", cc$lightSteelBlue(sp_surviv_region), "\n")
  
  return(hv_stats)
}