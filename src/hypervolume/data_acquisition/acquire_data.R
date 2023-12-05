source_all("./src/hypervolume/data_acquisition/components")

acquire_region = function(shapefiles, projection) {
  cat(blue("Acquiring regions. \n"))
  
  regions <- import_regions(shapefiles, projection, "./outputs/data_acquisition/region/logs/")
  
  return(regions)
}

reproject_region <- function(region, show_plot = F) {
  cat(blue("Reprojecting", strsplit(deparse(substitute(region)), "\\$")[[1]][[2]]), "\n")
  
  cat("Getting extents. \n")
  ext_east <- terra::ext(ext(region)$xmin, 0, ext(region)$ymin, ext(region)$ymax)
  ext_west <- terra::ext(0.00001, ext(region)$xmax, ext(region)$ymin, ext(region)$ymax)

  cat("Cropping in half. \n")
  vect_east <- terra::crop(region, ext_east)
  vect_west <- terra::crop(region, ext_west)

  cat("Reprojecting to longlat. \n")
  proj_east <- terra::project(vect_east, crs(longlat_crs))
  proj_west <- terra::project(vect_west, crs(longlat_crs))

  region_longlat <- rbind(proj_west, proj_east)
  
  if (show_plot == T) plot(region_longlat)
  
  cat(cc$lightGreen(strsplit(deparse(substitute(region)), "\\$")[[1]][[2]], "reprojected successfully. \n"))
  
  return(region_longlat)
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
