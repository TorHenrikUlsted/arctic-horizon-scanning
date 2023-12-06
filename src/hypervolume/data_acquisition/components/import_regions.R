import_regions <- function(shapefiles,
                           log_output,
                           verbose = T) {
  
  prj = crs("+proj=longlat")
  
  regions <- list()
  
  if (!dir.exists(log_output)) dir.create(log_output, recursive = T)
  
  log_output <- paste0(log_output, "regions_log.txt")
  
  appending <- F
  
  sink(log_output, append = appending)
  sink()
  
  cat(cc$aquamarine("Importing ", length(shapefiles), "shapefiles", "\n"))
  
  for (shapefile_name in names(shapefiles)) {
    if (is.na(shapefile_name) || shapefile_name == "") {
      cat(red("ERROR: You have to name your shapefiles, skipping missing a names... \n"))
      next
    }
    
    cat(cc$paleTurquoise("importing: ", shapefile_name, "\n"))
    
    shapefile <- shapefiles[shapefile_name]
    # Load the shapefiles
    region <- vect(shapefile)
    regions[[shapefile_name]] <- region
    
    
    # Get and print the extent
    ext_region <- ext(region)
    current_crs <- crs(region)
    
    if (verbose == T) {
      ext_region_printable <- as.vector(ext_region)
      cat("Extent of ", shapefile_name, ": ", ext_region_printable, "\n")
      
      # Check the original CRS
      original_crs <- crs(region)
      cat("Original CRS: ", crs(original_crs, proj = T), "\n")
    } else {
      next
    }
    
    # If the original CRS is not correctly defined, define it
    if (is.na(original_crs) || original_crs == "") {
      cat("Found blank or na crs, adding a placeholder crs. \n")
      crs(region) <- prj
    }
    
    # Append to log file if the file already exists
    if (file.exists(log_output)) appending = T else appending = F
    
    sink(log_output, append = appending)
    cat("Region: ", shapefile_name, "  |  ", "Extent: ", as.character(ext_region), "  |  ", "Projection: ", current_crs, "  |  ", "Append: ", appending, "\n")
    sink()
  }
  
  return(regions)
}