get_wc_data = function(var = "bio", res = 2.5, path = "./resources/region", version = "2.1", show.plot = F, verbose = F) {
  biovars <- list()
  catn("Getting World Clim data.")
  vebcat("Checking if files exist.", veb = verbose)
  # Check if all 19 biovariables are present, if not install
  for (i in 1:19) {
    dir_path <- paste0(path, "/wc", version, "_", res, "m/")
    file_name <- paste0(dir_path, "wc2.1_", res, "m_", var, "_", i, ".tif")
    zip_file <- paste0(dir_path, "wc2.1_", res, "m_", var, ".zip")
    
    if (!file.exists(file_name)) {
      vebcat(file_name, " is missing from the WorldClim data, looking for zip file...", color = "nonFatalError")
      
      if (!file.exists(zip_file)) {
        catn("Could not find zip file, trying to download...")
 
      } else {
        catn("Zip file located, unzipping...")
        unzip(zip_file, exdir = dir_path)
      }
      
      tryCatch(
        {
          r <- worldclim_global(var, res, path, version)
          
          if (!is.null(r)) {
            vebcat("WorldClim data downloaded successfully.", color = "proSuccess")
            biovars[[length(biovars) + 1]] <- r
          } else {
            vebcat("Failed to download worldClim data, is the server down?", color = "fatalError")
            message("Taking you to the site...")
            browseURL("https://www.worldclim.org/data/worldclim21.html")
            stop("Stopping, download failed")
          }
        },
        error = function(e) {
          vebcat("Failed to download worldClim data", color = "fatalError")
          stop("Stopping, download failed")
        }
      )
    } else {
      r <- rast(file_name)
      biovars[[length(biovars) + 1]] <- r
    }
  }
  
  if (show.plot == T) {
    for (i in 1:length(biovars)) {
      vebcat(paste0("Plotting bio_", as.character(i)), veb = verbose)
      # Plot each raster
      plot(biovars[[i]], main = paste("Biovar", i))
    }
  } else {
    vebcat("Plotting skipped.", veb = verbose)
  }
  
  biovars <- terra::rast(biovars)
  
  
  names(biovars) <- sub("wc2.1_2.5m_", "", (names(biovars)))
  varnames(biovars) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))
    
  return(biovars)
}