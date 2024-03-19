import_regions <- function(shapefiles, out.dir, verbose = F) {
  vebcat("Importing regions", color = "funInit")
  
  prj = crs(longlat_crs)
  
  regions <- list()
  
  create_dir_if(out.dir)
  
  log_output <- paste0(out.dir, "/regions-log.txt")
  
  appending <- F
  
  sink(log_output, append = appending)
  sink()
  
  if (is.character(shapefiles)) {
    shapefiles <- c(region = shapefiles)
  }
  
  vebcat("Importing ", length(shapefiles), "shapefile(s)", veb = verbose)
  
  for (shapefile_name in names(shapefiles)) {
    if (is.na(shapefile_name) || shapefile_name == "") {
      vebcat("ERROR: You have to name your shapefiles, skipping missing a names...", color = "nonFatalError")
      next
    }
    
    catn("importing: ", highcat(shapefile_name))
    
    shapefile <- shapefiles[shapefile_name]

    # Load the file using either rast or vect based on the file extension
    if (grepl("\\.shp$", shapefile)) {
      region <- terra::vect(shapefile)
    } else if (grepl("\\.tif$", shapefile) || grepl("\\.nc$", shapefile)) {
      region <- terra::rast(shapefile)
    } else {
      stop("Unsupported file type", shapefile)
    }

    # Get and print the extent
    ext_region <- ext(region)
    current_crs <- crs(region, proj = T)
    
    ext_region_printable <- as.vector(ext_region)
    vebcat("Extent of ", shapefile_name, ": ", ext_region_printable, veb = verbose)
    
    # Check the original CRS
    original_crs <- crs(region)
    vebcat("Original CRS: ", crs(original_crs, proj = T), veb = verbose)
    
    # If the original CRS is not correctly defined, define it
    if (is.na(original_crs) || original_crs == "") {
      vebcat("Found blank or na crs. Needs manual processing", color = "nonFatalError")
    }

    regions[[shapefile_name]] <- region
    
    # Append to log file if the file already exists
    if (file.exists(log_output)) appending = T else appending = F
    
    sink(log_output, append = appending)
    catn("Region: ", shapefile_name, "  |  ", "Extent: ", ext_region, "  |  ", "Projection: ", current_crs)
    sink()
  }
  
  vebcat("Regions imported successfully", color = "funSuccess")
  
  return(regions)
}
