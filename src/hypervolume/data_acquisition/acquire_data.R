source_all("./src/hypervolume/data_acquisition/components")

acquire_data <- function(shapefiles, projection, sp_df = NULL, plot_it, verbose) {
  cat(blue("Acquiring data. \n"))
  
  cat("importing regions. \n")
  if (file.exists("./outputs/data_acquisition/region/cavm/biovars_cavm.rds")) {
    cat("File found", green(O_O), ", Loading file. \n")
    
    wc_cavm <- readRDS("./outputs/data_acquisition/region/cavm/biovars_cavm.rds")
    
  } else {
    regions <- import_regions(shapefiles, projection, "./outputs/data_acquisition/region/logs/")
    
    region_crop <- start_timer("crop_region")
    cat("Cropping worldClim data to CAVM. \n")
    
    wc_cavm <- wc_to_region(regions$cavm, projection = projection, plot_it = plot_it, verbose = verbose)
    
    end_timer(region_crop)
    
    cat("Saving raster stack to:", cc$lightSteelBlue("./outputs/data_acquisition/region/cavm/biovars_cavm.rds"), "\n")
    
    saveRDS(wc_cavm, file = "./outputs/data_acquisition/region/cavm/biovars_cavm.rds")
    
  }
  
  # Leaving blank gets worldclim data for the whole world
  cat("Acquiring WorldClim data for the entire world. \n")
  wc_world = wc_to_region()
  
  cat("Acquiring species. \n")
  
  if (is.null(sp_df)) {
    cat(cc$lightCoral("Species data frame not found. \n"))
    
    sp_df <- setup_env(test = T)
  } else {
    cat("Species data frame found. \n")
    
  }
  
  cat(cc$lightGreen("Data Acquisition successfully completed. \n"))
  
  return(list(
    wc_cavm,
    wc_world,
    sp_df
  ))
}