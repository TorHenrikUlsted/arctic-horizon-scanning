source_all("./src/hypervolume/data_acquisition/components")

acquire_region = function(shapefiles) {
  cat(blue("Acquiring regions. \n"))
  
  regions <- import_regions(shapefiles, "./outputs/data_acquisition/region/logs/")
  
  return(regions)
}

acquire_biovars = function(show_plot = F) {
  cat(blue("Acquiring WorldClim biovariables. \n"))
  
  biovars <- get_wc_data(show_plot = show_plot)
  
  cat(cc$lightGreen("Biovars acquired successfully \n"))
  return(biovars)
}

acquire_region_data = function(biovars, regions, projection, show_plot = F, verbose = T) {
  cat(blue("Acquiring WorldClim region data. \n"))
  
  region_data_list <- list()
  
  for(i in seq_along(regions)) {
    region <- regions[[i]]
    
    name <- names(regions)[i]
    
    cat("Using", cc$lightSteelBlue(name), "\n")
    
    filename <- paste0("./outputs/data_acquisition/region/", name, "/biovars_", projection, ".rds") 
    
    
    if (file.exists(filename)) {
    cat("File found", green("O_O"), "Loading file... \n")
    
    region_data <- readRDS(filename)
    
    } else {
      cat(red("File not found. \n"))
      
      region_crop <- start_timer("crop_region")
      
      region_data <- wc_to_region(biovars, region, projection, show_plot = F, verbose = T)
      
      end_timer(region_crop)
      
      cat("Saving", cc$lightSteelBlue(name), "to:", cc$lightSteelBlue(filename), "\n")
      
      create_dir_if(paste0("./outputs/data_acquisition/region/", name))
      
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

acquire_species_data = function(sp_df, test, big_test) {
  cat("Acquiring species. \n")
  
  if (is.null(sp_df)) {
    cat(cc$lightCoral("Species data frame not found. \n"))
    
    sp_df <- setup_sp(test = test, big_test = big_test)
    
  } else cat("Species data frame found. \n")

  cat(cc$lightGreen("Species data acquired successfully. \n"))
  
  return(sp_df)
}
