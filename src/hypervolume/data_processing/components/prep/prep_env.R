prepare_environment <- function(sp_points, biovars, verbose = T) {
  
  # Create an empty matrix to store environmental values
  env_values <- matrix(nrow = nrow(sp_points), ncol = terra::nlyr(biovars))
  colnames(env_values) <- names(biovars)
  #set the crs()
  crs(sp_points) <- crs(biovars[[1]])
  
  if (verbose == T) {
    if (identical(terra::crs(sp_points), terra::crs(biovars[[1]]))) {
      cat("CRS for bio variables and species are identical:", green("TRUE"), "\n") 
    } else cat("CRS for bio variables and species are identical:", red("FALSE"), "\n")
  }
  
  # Extract environmental values for each point
  for (i in 1:terra::nlyr(biovars)) {
    if (verbose == T) {
      cat("Extracting biovariable for:", cc$lightSteelBlue(names(biovars)), "\n")
    }
    
    env_values[, i] <- terra::extract(biovars[[i]], sp_points, df=TRUE, ID = FALSE)[,1]
  }
  

  if (verbose) cat("Cleaning extracted data. \n")
  
  if (any(is.na(env_values)) || any(env_values == "")) {
    if (verbose) cat("Some extracted values are NA:", red("TRUE"), "\n")
    
    env_values[env_values == ""] <- NA
    
    clean_env_matrix <- na.omit(env_values)
    
    if (verbose == T) {
      if (any(is.na(clean_env_matrix))) cat(red("Cleaning NA-values failed. \n")) else cat(green("Cleaning NA-values was successfull. \n")) 
    }

    
  }  else {
    if (verbose) cat("Some extracted values are NA:", green("FALSE"), "\n")
    
    clean_env_matrix <- env_values
  } 
  
  return(clean_env_matrix)
}

