get_region_data = function(biovars, regions, projection, show_plot = F, verbose = F) {
  vebcat("Acquiring WorldClim region data.", color = "funInit")
  
  region_data_list <- list()
  
  for(i in seq_along(regions)) {
    region <- regions[[i]]
    
    name <- names(regions)[i]
    
    catn("Using", highcat(name))
    
    dir_name <- paste0("./outputs/hypervolume/data_acquisition/region/", name)
    create_dir_if(dir_name)
    
    filename <- paste0(dir_name, "/biovars_", projection, ".rds") 
    
    
    if (file.exists(filename)) {
      vebcat("File found", "Loading file...", veb = verbose)
      
      region_data <- readRDS(filename)
      
    } else {
      vebcat("File not found.", color = "nonFatalError", veb = verbose)
      
      region_crop <- start_timer("crop_region")
      
      region_data <- wc_to_region(biovars, region, projection, show_plot = F, verbose)
      
      end_timer(region_crop)
      
      catn("Saving", highcat(name), "to:", colcat(filename, color = "output"))
      
      saveRDS(region_data, file = filename)
    }
    # Add the region data to the list
    region_data_list[[name]] <- region_data
  }
  
  vebcat("Region data acquired successfully.", color = "funSuccess")
  
  # Check if the list has only one item
  if(length(region_data_list) == 1) {
    return(region_data_list[[1]])
  } else {
    return(region_data_list)
  }
}