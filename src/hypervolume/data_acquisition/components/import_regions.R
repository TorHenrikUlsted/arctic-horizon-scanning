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
    
    
    if (projection == "longlat") {
      
      cat("Choosing", cc$lightSteelBlue("longlat"),  "coordinate system. \n")
      if (grepl("+proj=longlat", original_crs, fixed = TRUE) == FALSE) cat(red("Is not long Lat, needs reprojection. \n"))
      print(crs(longlat_crs, proj = T))
      prj <- longlat_crs
      
    } else if (projection == "laea") {
      
      cat("Choosing", cc$lightSteelBlue("polar"), "coordinate system. \n")
      if (grepl("+proj=laea", crs(original_crs, proj = T), fixed = TRUE) == FALSE) cat(red("Is not polar \n"))
      prj <- cavm_crs
      
    } else stop("You can only choose projection 'longlat' or 'laea'.")
    

    if (!isTRUE(identical(original_crs, prj))) {
      if (verbose == T) cat(yellow("Original CRS not identical to current CSR. \n"))
      if (verbose == T) cat("Reprojecting", cc$lightSteelBlue(shapefile_name), "to: ", crs(prj, proj = T), "\n") else next
      
      regions[[shapefile_name]] <- terra::project(region, prj)
      
      if (!isTRUE(identical(current_crs, prj))) cat(green("Reprojection completed successfully: "), crs(regions[[shapefile_name]], proj = T), "\n") else cat(red("Reprojection failed."))
    } else {
      if (verbose == T) cat("Original CRS identical to current CSR. \n")
    }
    
    # Append to log file if the file already exists
    if (file.exists(log_output)) appending = T else appending = F
    
    sink(log_output, append = appending)
    cat("Region: ", shapefile_name, "  |  ", "Extent: ", as.character(ext_region), "  |  ", "Projection: ", current_crs, "  |  ", "Append: ", appending, "\n")
    sink()
  }
  
  return(regions)
}