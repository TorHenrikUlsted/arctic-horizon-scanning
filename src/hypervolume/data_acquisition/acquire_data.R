source_all("./src/hypervolume/data_acquisition/components")

acquire_region_data = function(biovars, regions, projection, show_plot = F, verbose = F) {
  cat(blue("Acquiring WorldClim region data. \n"))
  
  region_data_list <- list()
  
  for(i in seq_along(regions)) {
    region <- regions[[i]]
    
    name <- names(regions)[i]
    
    cat("Using", cc$lightSteelBlue(name), "\n")
    
    dir_name <- paste0("./outputs/hypervolume/data_acquisition/region/", name)
    create_dir_if(dir_name)
    
    filename <- paste0(dir_name, "/biovars_", projection, ".rds") 
    
    
    if (file.exists(filename)) {
      if (verbose) cat("File found", "Loading file... \n")
    
    region_data <- readRDS(filename)
    
    } else {
      if (verbose) cat(red("File not found. \n"))
      
      region_crop <- start_timer("crop_region")
      
      region_data <- wc_to_region(biovars, region, projection, show_plot = F, verbose)
      
      end_timer(region_crop)
      
      if (verbose) cat("Saving", cc$lightSteelBlue(name), "to:", cc$lightSteelBlue(filename), "\n")
      
      saveRDS(region_data, file = filename)
    }
    # Add the region data to the list
    region_data_list[[name]] <- region_data
  }
  
  cat(cc$lightGreen("Region data acquired successfully. \n"))
  
  # Check if the list has only one item
  if(length(region_data_list) == 1) {
    return(region_data_list[[1]])
  } else {
    return(region_data_list)
  }
}
