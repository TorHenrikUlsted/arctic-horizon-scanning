import_regions <- function(shapefiles,
                           log_output,
                           verbose = F) {
  log_output = paste0(log_output, "/")
  
  prj = crs(longlat_crs)
  
  regions <- list()
  
  if (!dir.exists(log_output)) dir.create(log_output, recursive = T)
  
  log_output <- paste0(log_output, "regions_log.txt")
  
  appending <- F
  
  sink(log_output, append = appending)
  sink()
  
  if (verbose) cat(cc$aquamarine("Importing ", length(shapefiles), "shapefiles", "\n"))
  
  for (shapefile_name in names(shapefiles)) {
    if (is.na(shapefile_name) || shapefile_name == "") {
      cat(red("ERROR: You have to name your shapefiles, skipping missing a names... \n"))
      next
    }
    
    cat(cc$paleTurquoise("importing: ", shapefile_name, "\n"))
    
    shapefile <- shapefiles[shapefile_name]
    # Get the file extensions

    # Load the shapefiles
    # Load the file using either rast or vect based on the file extension
    if (grepl("\\.shp$", shapefile)) {
      region <- vect(shapefile)
    } else if (grepl("\\.tif$", shapefile) || grepl("\\.nc$", shapefile)) {
      region <- rast(shapefile)
    } else {
      stop("Unsupported file type", shapefile)
    }
    
    # Get and print the extent
    ext_region <- ext(region)
    current_crs <- crs(region)
    
    ext_region_printable <- as.vector(ext_region)
    if (verbose) cat("Extent of ", shapefile_name, ": ", ext_region_printable, "\n")
    
    # Check the original CRS
    original_crs <- crs(region)
    if (verbose) cat("Original CRS: ", crs(original_crs, proj = T), "\n")

    
    # If the original CRS is not correctly defined, define it
    if (is.na(original_crs) || original_crs == "") {
      cat(red("Found blank or na crs. Needs manual processing. \n"))
    }
  
    regions[[shapefile_name]] <- region
    
    # Append to log file if the file already exists
    if (file.exists(log_output)) appending = T else appending = F
    
    sink(log_output, append = appending)
    cat("Region: ", shapefile_name, "  |  ", "Extent: ", as.character(ext_region), "  |  ", "Projection: ", current_crs, "  |  ", "Append: ", appending, "\n")
    sink()
  }
  
  return(regions)
}