analyze_hv_stats <- function(region_hv, sp_hv, spec.name, verbose) {
  catn("Analyzing hypervolume for", highcat(sp_hv@Name), "species. \n")
  
  hv_set <- hypervolume_set(sp_hv, region_hv, check.memory = F)
  
  hv_stats <- hypervolume_overlap_statistics(hv_set)
  
  sp_surviv_region <- 1 - hv_stats[[4]]
  
  catn("Volume of", spec.name, "overlap in the CAVM", highcat(sp_surviv_region))
  
  return(hv_stats)
}

analyze_region_hv <- function(biovars, out.dir, name, method, verbose) {
  vebcat("Analyzing cavm hypervolume.", color = "funInit")
  
  directory <- paste0(out.dir, "/", tolower(name), "/")
  
  create_dir_if(directory)
  
  if (!file.exists(paste0(directory, "hypervolume_", method, ".rds"))) {
    vebcat("File not found, initating hypervolume sequence.", veb = verbose)
    matrix <- terra::values(biovars, na.rm = T)
    
    if (any(is.na(matrix))) {
      vebcat("Some biovars region values are NA.", color = "nonFatalError")
    } else {
      vebcat("No biovars region values are NA.", color = "proSuccess")
    }
    
    vebprint(head(matrix, 3), text = "Matrix sample:")
    
    hv <- hypervolume(
      matrix, 
      name = name,
      method = method, 
      verbose = verbose
    )
    
    hv_save <- paste0(directory, "hypervolume_", method, ".rds")
    
    catn("Saving region hypervolume to:", colcat(hv_save, color = "output"))    
    
    saveRDS(hv, hv_save)
  } else {
    catn(name, "Hypervolume found, loading file.")
    # Load the hypervolume from the file
    hv <- readRDS(paste0(directory, "hypervolume_", method, ".rds"))
  }
  
  vebcat(name, paste0("Hypervolume_", method), "Analysis Complete.", color = "funSuccess")
  
  return(hv)
}

analyze_ind_hv <- function(env_data, verbose = FALSE) {
  vebcat("Analyzing hypervolumes for", highcat(length(env_data)), "species.", color = "funInit")
  
  hv_list <- list()
  
  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(env_data), style = 3)
  
  for (species in names(env_data)) {
    catn("Creating hypervolume for species number", highcat(species))
    # Get the data for the species
    data <- env_data[[species]]
    
    vbprint(head(data, 3), veb = verbose)
    # Initialize an empty list to store the hypervolumes for the species
    species_hvs <- list()
    
    # For each dimension in the data
    for (dimension in colnames(data)) {
      catn("Using dimension:", highcat(dimension))
      print(head(data[, dimension]), 3)
      # Create a hypervolume for the dimension
      hv <- hypervolume_box(data[, dimension, drop = FALSE])
      catn("Add", dimension, "to species list.")
      # Add the hypervolume to the list
      species_hvs[[dimension]] <- hv
      
      vebcat("Adding dimension to list.", veb = verbose)
    }
    
    # Create a HypervolumeList for the species
    catn("Creating hypervolumeList")
    species_hv_list <- do.call(hypervolume_join, species_hvs)
    
    # Add the HypervolumeList to the list
    vebcat("Add species hypervolumeList to main list.", veb = verbose)
    hv_list[[species]] <- species_hv_list
    
    # Update the progress bar
    setTxtProgressBar(pb, which(names(env_data) == species))
  }; catn()
  
  # Close the progress bar
  close(pb)
  
  vebcat("Analyzed hypervolumes for", highcat(length(env_data)), "species successfully", color = "funSuccess")
  
  return(hv_list)
}