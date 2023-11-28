source_all("./src/hypervolume/data_processing/components")

process_data <- function(biovars, sp_data, projection, verbose) {
  cat(blue("Initiating data processing protocol. \n"))
  
  # choose wanted variables
  cat("Using biovars:", cc$lightSteelBlue(sub("wc2.1_2.5m_", "", as.character(names(biovars)))), "\n")
  
  sp_list <- unique(sp_data$species)
 
  cat("Preparing", cc$lightSteelBlue(length(sp_list)), "species. \n")
  
  sp_points_list <- list()
  
  for (i in seq_along(sp_list)) {
    cat("\nPreparing", cc$lightSteelBlue(sp_list[i]), "\n")
    # Subset the data for the current species
    df_species <- sp_data[sp_data$species == sp_list[i], ]
    
    # Prepare the species occurrence data for the current species
    sp_points <- prepare_species(df_species, projection, verbose)
    
    # Store the spatial points data frame in the list
    sp_points_list[[sp_list[i]]] <- sp_points
  }

  cat("Preparing", cc$lightSteelBlue(length(sp_points_list)), "environments. \n")
  
  env_values_list <- list()
  
  # Loop over each unique species
  for (i in seq_along(sp_points_list)) {
    # Get the spatial points data frame for the current species
    sp_points_species <- sp_points_list[[i]]
    
    # Prepare the environmental data for the current species
    env_values_species <- prepare_environment(sp_points_species, biovars, verbose)
    
    # Store the environmental data in the list
    env_values_list[[names(sp_points_list)[i]]] <- env_values_species
  }
  
  #env_values_list <- check_zero_var(env_values_list, verbose = verbose)
  
  cat(cc$lightGreen("data processing sucessfully completed. \n"))

  return(env_values_list)
}
