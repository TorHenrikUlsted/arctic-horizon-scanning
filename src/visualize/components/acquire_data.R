acquire_visualize_data <- function(hv.dir) {
  # Get species stats

  
  # Get inclusion tif
  sp_dirs <- list.dirs("./outputs/hypervolume/sequence/projections")
  
  inc_stack <- terra::rast()
  
  for (i in sp_dirs) {
    sp_inclusion <- paste0(sp_dirs[[i]], "/inclusion-0.5.tif")
    
    sp_rast <- terra::rast(sp_inclusion)
    
    inc_stack <- c(inc_stack, sp_rast)
  }
  
  
  # Get probability tif
  
  prob_stack <- terra::rast()
  
  for (i in sp_dirs) {
    sp_prob <- paste0(sp_dirs[[i]], "/probability.tif")
    
    sp_rast <- terra::rast(sp_prob)
    
    prob_stack <- c(prob_stack, sp_rast)
  }
  
  
  
  
  return(list(
    included_sp = included_sp,
    excluded_sp = excluded_sp,
    inc_stack = inc_stack,
    prob_stack = prob_stack
  ))
}