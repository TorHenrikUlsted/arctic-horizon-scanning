source_all("./src/hypervolume/data_processing/components")

scale_biovars <- function(biovars, verbose = F) {
 if (verbose) cat(blue("Scaling biovariables \n"))
  name = deparse(substitute(biovars))
  
  save_dir <- "./outputs/hypervolume/data_processing/region"
  create_dir_if(save_dir)
  
  save_path <- paste0(save_dir, "/", name, "_scaled.rds") 
  
  # Check if the scaled data already exists
  if (file.exists(save_path)) {
    if (verbose) cat("Scaled data found, loading from file. \n")
    scaled_biovars <- readRDS(save_path)
    if (verbose) cat(cc$lightGreen("Scaling biovariables completed successfully. \n"))
    return(scaled_biovars)
  }
  
  scaled_biovars <- biovars
  # Scale the dimensions
  for (i in 1:nlyr(scaled_biovars)) {
    cat("Scaling", cc$lightSteelBlue(names(scaled_biovars)[i]), "\n")
    # Extract the i-th layer
    layer <- biovars[[i]]
    
    # Scale the layer
    layer_scaled <- app(layer, fun = scale)
    
    # Replace the i-th layer of r_scaled with the scaled layer
    scaled_biovars[[i]] <- layer_scaled
  }
  
  if (verbose) cat("Renaming layers. \n")
  names(scaled_biovars) <- names(biovars)
  
  # Save the scaled data
  saveRDS(scaled_biovars, save_path)
  cat(name, "scaled data saved to file. \n")
  
  if (verbose) cat(cc$lightGreen("Scaling biovariables completed successfully. \n"))
  return(scaled_biovars)
}
