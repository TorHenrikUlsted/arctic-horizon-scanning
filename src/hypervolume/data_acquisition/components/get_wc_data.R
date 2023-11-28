get_wc_data = function(region = NULL, var = "bio", res = 2.5, path = "./resources/region", version = "2.1") {
  biovars <- list()
  cat("Checking if files exist \n")
  # Check if all 19 biovariables are present, if not install
  for (i in 1:19) {
    dir_path <- paste0(path, "/wc", version, "_", res, "m/")
    file_name <- paste0(dir_path, "wc2.1_", res, "m_", var, "_", i, ".tif")
    zip_file <- paste0(dir_path, "wc2.1_", res, "m_", var, ".zip")
    
    if (!file.exists(file_name)) {
      cat(red(file_name, " is missing from the WorldClim data, looking for zip file... \n"))
      
      if (!file.exists(zip_file)) {
        cat("Could not find zip file, trying to download... \n")
        next
      } else {
        cat("Zip file located O_O, unzipping... \n")
        unzip(zip_file, exdir = dir_path)
        
        next
      }
      
      tryCatch(
        {
          r <- worldclim_global(var, res, path, version)
          
          if (!is.null(r)) {
            cat(green("WorldClim data downloaded successfully. \n"))
          } else {
            cat(red("Failed to download worldClim data, is the server down? \n"))
            message("Taking you to the site...")
            browseURL("https://www.worldclim.org/data/worldclim21.html")
            stop("Exiting, download failed")
          }
        },
        error = function(e) {
          cat(red("Failed to download worldClim data \n"))
          stop("Exiting, download failed")
        }
      )
    } else {
      r <- rast(file_name)
    }
    
    biovars[[length(biovars) + 1]] <- r
  }
  
  return(biovars)
}