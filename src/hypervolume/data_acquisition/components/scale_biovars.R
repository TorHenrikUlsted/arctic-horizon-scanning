scale_biovars <- function(biovars, verbose = F) {
  name = deparse(substitute(biovars))
  vebcat("Scaling biovariables", name, color = "funInit")
  
  save_dir <- "./outputs/hypervolume/data_processing/region"
  create_dir_if(save_dir)
  
  save_path <- paste0(save_dir, "/", name, "_scaled.rds")
  vebcat("save_path:", save_path, veb = verbose)
  
  # Check if the scaled data already exists
  if (file.exists(save_path)) {
    catn("Scaled data found, loading from file.")
    scaled_biovars <- readRDS(save_path)
    vebcat("Scaling biovariables completed successfully.", color = "funSuccess")
    return(scaled_biovars)
  }
  
  scaled_biovars <- biovars
  # Scale the dimensions
  for (i in 1:nlyr(scaled_biovars)) {
    catn("Scaling", highcat(names(scaled_biovars)[i]))
    # Extract the i-th layer
    layer <- biovars[[i]]
    
    # Scale the layer
    layer_scaled <- app(layer, fun = scale)
    
    # Replace the i-th layer of r_scaled with the scaled layer
    scaled_biovars[[i]] <- layer_scaled
  }
  
  vebcat("Renaming layers.", veb = verbose)
  names(scaled_biovars) <- names(biovars)
  
  # Save the scaled data
  saveRDS(scaled_biovars, save_path)
  catn(name, "scaled data saved to file:", colcat(save_path, color = "output"))
  
  vebcat("Scaling biovariables completed successfully.", color = "funSuccess")
  return(scaled_biovars)
}