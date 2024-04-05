get_region_data = function(biovars, shapefile, projection, show_plot = F, verbose = F) {
  vebcat("Acquiring WorldClim region data.", color = "funInit")
  
  region <- load_region(shapefile, verbose = verbose)
    
  name <- basename(shapefile)
  name <- head(strsplit(name, split = "\\.")[[1]], -1)
  
  catn("Using", highcat(name))
  
  dir_name <- paste0("./outputs/hypervolume/data_acquisition/region/", name)
  create_dir_if(dir_name)
  
  filename <- paste0(dir_name, "/biovars_", projection, ".rds") 
  
  
  if (file.exists(filename)) {
    vebcat("File found", "Loading file...", veb = verbose)
    
    region_data <- readRDS(filename)
    
  } else {
    vebcat("File not found.", color = "nonFatalError", veb = verbose)
    
    region_crop_timer <- start_timer("crop_region")
    
    region_data <- wc_to_region(biovars, region, projection, show_plot = F, verbose)
    
    end_timer(region_crop_timer)
    
    catn("Saving", highcat(name), "to:", colcat(filename, color = "output"))
    
    saveRDS(region_data, file = filename)
  }


  vebcat("Region data acquired successfully.", color = "funSuccess")

# Check if the list has only one item
  return(region_data)
}
