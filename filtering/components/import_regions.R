import_regions = function(shapefiles, projection, log_output = "./outputs/regions.txt") {
  regions = list()
  appending = F
  
  sink(log_output, append = appending)
  sink()
  
  for (shapefile_name in names(shapefiles)) {
    cat("Importing shapefile: ", shapefile, "\n")
    
    shapefile = shapefiles[shapefile_name]
    
    # Load the shapefiles
    region = vect(shapefile)
    
    # Get and print the extent
    ext_region = ext(region)
    cat("Extent of ", shapefile, ": ", ext_region, "\n")
    
    current_crs = crs(region, proj = T, describe = T)
    cat("Current CRS: ", current_crs, "\n")
    
    if (current_crs != projection) {
      cat("Reprojecting ", shapefile, " to: ", projection, "\n")
      region = terra::project(region, projection)
    }
    
    if (!is.null(shapefile_name)) region_name = shapefile_name else region_name = basename(shapefile)
    regions[[region_name]] = region
    
    # Append to log file if the file already exists
    if (file.exists(log_output)) appending = T else appending = F
    
    sink(log_output, append = appending)
    cat("Region: ", region, "  |  ", "Extent: ", ext_region, "  |  ", "Append: ", appending)
    sink()
  }
  
  return(regions)
}