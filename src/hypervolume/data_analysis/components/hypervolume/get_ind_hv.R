get_ind_hv <- function(env_data) {
  cat(blue("Analyzing hypervolumes for"), cc$lightSteelBlue(length(env_data)), blue("species. \n"))
  
  hv_list <- list()
  
  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(env_data), style = 3)
  
  for (species in names(env_data)) {
    cat("Creating hypervolume for species number", cc$lightSteelBlue(species), "\n")
    # Get the data for the species
    data <- env_data[[species]]
    
    print(head(data, 3))
    # Initialize an empty list to store the hypervolumes for the species
    species_hvs <- list()
    
    # For each dimension in the data
    for (dimension in colnames(data)) {
      cat("Using dimension:", cc$lightSteelBlue(dimension), "\n")
      print(head(data[, dimension]), 3)
      # Create a hypervolume for the dimension
      hv <- hypervolume_box(data[, dimension, drop = FALSE])
      cat("Add", dimension, "to species list. \n")
      # Add the hypervolume to the list
      species_hvs[[dimension]] <- hv
      
      cat("Adding dimension to list \n")
    }
    
    # Create a HypervolumeList for the species
    cat("Creating hypervolumeList \n")
    species_hv_list <- do.call(hypervolume_join, species_hvs)
    
    # Add the HypervolumeList to the list
    cat("Add species hypervolumeList to main list \n")
    hv_list[[species]] <- species_hv_list
    
    # Update the progress bar
    setTxtProgressBar(pb, which(names(env_data) == species))
    cat("\n")
  }
  
  # Close the progress bar
  close(pb)
  
  return(hv_list)
}