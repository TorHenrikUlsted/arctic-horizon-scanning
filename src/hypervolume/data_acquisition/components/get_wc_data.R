get_wc_data = function(var = "bio", res = 2.5, path = "./resources/region", version = "2.1", show.plot = F, verbose = F) {
  biovars <- list()
  if (verbose) cat("Checking if files exist \n")
  # Check if all 19 biovariables are present, if not install
  for (i in 1:19) {
    dir_path <- paste0(path, "/wc", version, "_", res, "m/")
    file_name <- paste0(dir_path, "wc2.1_", res, "m_", var, "_", i, ".tif")
    zip_file <- paste0(dir_path, "wc2.1_", res, "m_", var, ".zip")
    
    if (!file.exists(file_name)) {
      cat(red(file_name, " is missing from the WorldClim data, looking for zip file... \n"))
      
      if (!file.exists(zip_file)) {
        cat("Could not find zip file, trying to download... \n")
 
      } else {
        cat("Zip file located, unzipping... \n")
        unzip(zip_file, exdir = dir_path)
      }
      
      tryCatch(
        {
          r <- worldclim_global(var, res, path, version)
          
          if (!is.null(r)) {
            cat(green("WorldClim data downloaded successfully. \n"))
            biovars[[length(biovars) + 1]] <- r
          } else {
            cat(red("Failed to download worldClim data, is the server down? \n"))
            message("Taking you to the site...")
            browseURL("https://www.worldclim.org/data/worldclim21.html")
            stop("Stopping, download failed")
          }
        },
        error = function(e) {
          cat(red("Failed to download worldClim data \n"))
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
      if (verbose) cat(paste0("Plotting bio_", as.character(i)), "\n")
      # Plot each raster
      plot(biovars[[i]], main = paste("Biovar", i))
    }
  } else {
    if (verbose) cat("Plotting skipped \n")
  }
  
  biovars <- terra::rast(biovars)
  
  
  names(biovars) <- sub("wc2.1_2.5m_", "", (names(biovars)))
  varnames(biovars) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))
    
  return(biovars)
}