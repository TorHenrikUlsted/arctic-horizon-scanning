source_all("./src/hypervolume/data_processing/components")

scale_biovars <- function(biovars) {
  cat(blue("Scaling biovariables \n"))
  
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
  
  cat("Renaming layers. \n")
  names(scaled_biovars) <- names(biovars)
  
  cat(cc$lightGreen("Scaling biovariables completed successfully. \n"))
  return(scaled_biovars)
}

process_sp_data <- function(sp_data, projection, verbose) {
  cat(blue("Processing species data.\n"))
  
  sp_list <- unique(sp_data$species)
  
  sp_points_list <- list()
  
  cat("Preparing", cc$lightSteelBlue(length(sp_list)), "species. \n")
  
  for (i in seq_along(sp_list)) {
    cat("\nPreparing", cc$lightSteelBlue(sp_list[i]), "\n")
    # Subset the data for the current species
    df_species <- sp_data[sp_data$species == sp_list[i], ]
    
    # Prepare the species occurrence data for the current species
    sp_points <- prepare_species(df_species, projection, verbose)
    
    # Store the spatial points data frame in the list
    sp_points_list[[sp_list[i]]] <- sp_points
  }
  
  cat(cc$lightGreen("Species preperation competed successfully. \n"))
  
  return(sp_points_list)
}

process_env_data <- function(biovars, sp_data, projection, verbose) {
  cat(blue("Processing environment data. \n"))
  
  cat("Using biovars:", cc$lightSteelBlue(sub("wc2.1_2.5m_", "", as.character(names(biovars)))), "\n")

  cat("Preparing", cc$lightSteelBlue(length(sp_data)), "environments. \n")
  
  env_values_list <- list()
  
  # Loop over each unique species
  for (i in seq_along(sp_data)) {
    # Get the spatial points data frame for the current species
    sp_points_species <- sp_data[[i]]
    
    # Prepare the environmental data for the current species
    env_values_species <- prepare_environment(sp_points_species, biovars, verbose)
    
    # Store the environmental data in the list
    env_values_list[[names(sp_data)[i]]] <- env_values_species
    
  }
  
  cat(cc$lightGreen("Environment preperation competed successfully. \n"))

  return(env_values_list)
}
