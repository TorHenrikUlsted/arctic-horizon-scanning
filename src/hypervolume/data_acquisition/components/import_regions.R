import_regions <- function(shapefiles,
                           projection = "longlat",
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
    current_crs <- crs(region, proj = T)
    
    if (verbose == T) {
      ext_region_printable <- as.vector(ext_region)
      cat("Extent of ", shapefile_name, ": ", ext_region_printable, "\n")
      
      # Check the original CRS
      original_crs <- crs(region, proj = T)
      cat("Original CRS: ", original_crs, "\n")
    } else {
      next
    }
    
    
    # If the original CRS is not correctly defined, define it
    if (is.na(original_crs) || original_crs == "") {
      cat("Found blank or na crs, adding a placeholder crs. \n")
      crs(region) <- prj
    }
    
    if (projection == "longlat") {
      cat("Choosing", cc$lightSteelBlue("longlat"),  "coordinate system. \n")
      
      if (grepl("+proj=longlat", original_crs, fixed = TRUE) == FALSE) {
        cat(red("Is not long Lat \n"))
        prj <- crs("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84")
      } else {
        prj <- crs("+proj=longlat +datum=WGS84")
      }
    } else if (projection == "laea") {
      cat("Choosing", cc$lightSteelBlue("polar"), "coordinate system. \n")
      
      if (grepl("+proj=laea", original_crs, fixed = TRUE) == FALSE) {
        cat(red("Is not polar \n"))
        prj <- crs("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
      } else {
        prj <- crs("+proj=laea +datum=WGS84")
      }
    } else {
      stop("You can only choose projection 'longlat' or 'polar'.")
    }
    

    if (!isTRUE(identical(original_crs, prj))) {
      if (verbose == T) cat("Reprojecting ", shapefile_name, " to: ", crs(prj, proj = T), "\n") else next
      
      regions[[shapefile_name]] <- terra::project(region, prj)
      
      if (!isTRUE(identical(current_crs, prj))) cat(green("Reprojection completed successfully: "), crs(regions[[shapefile_name]], proj = T), "\n") else cat(red("Reprojection failed."))
    } else {
      cat("Reprojection not needed")
    }
    
    
    # Append to log file if the file already exists
    if (file.exists(log_output)) appending = T else appending = F
    
    sink(log_output, append = appending)
    cat("Region: ", shapefile_name, "  |  ", "Extent: ", as.character(ext_region), "  |  ", "Projection: ", current_crs, "  |  ", "Append: ", appending, "\n")
    sink()
  }
  
  return(regions)
}