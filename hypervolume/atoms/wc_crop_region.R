wc_to_region = function(region, var = "bio", res = 2.5, path = "./resources", version = "2.1", plot = T, verbose= T) {
  biovars = list()
  cat("Checking if files exist \n")
  # Check if all 19 biovariables are present, if not install
  for (i in 1:19) {
    dir_path = paste0(path,"/wc",version,"_",res,"m/")
    file_name = paste0(dir_path,"wc2.1_",res,"m_",var,"_",i,".tif")
    zip_file = paste0(dir_path,"wc2.1_",res,"m_",var,".zip")
    
    if(!file.exists(file_name)) {
      cat(red(file_name, " is missing from the WorldClim data, looking for zip file... \n"))
      
      if(!file.exists(zip_file)) {
        cat("Could not find zip file, trying to download... \n")
        next
      } else {
        cat("Zip file located O_O, unzipping... \n")
        unzip(zip_file, exdir = dir_path)
        
        next
      }
      
      tryCatch({
        r = worldclim_global(var, res, path, version)
        
        if (!is.null(r)) {
          cat(green("WorldClim data downloaded successfully. \n"))
        } else {
          cat(red("Failed to download worldClim data, is the server down? \n"))
          message("Taking you to the site...")
          browseURL("https://www.worldclim.org/data/worldclim21.html")
          stop("Exiting, download failed")
        }
      }, error = function(e) {
        cat(red("Failed to download worldClim data \n"))
        stop("Exiting, download failed")
      })
      
    } else {
      r = rast(file_name)
    }
    
    biovars[[length(biovars) + 1]] = r
    
  }
  
  cat("Starting cropping sequence \n")
  
  biovarsMask = list()
  
  for (biovar in biovars) {
    if (!isTRUE(identical(crs(biovar), crs(region) ))) {
      if (verbose == T){
        cat(red("Crs not identical: \n"),
            "Biovars: ", as.character(crs(biovar, proj = T, describe = T)), "\n",
            "Region", as.character(crs(region, proj = T, describe = T)), "\n",
            "Making them identical...", "\n"
        )
        
        crs(region) = crs(biovar)
        
        cat(green("Crs fixed: \n"),
            "Biovars: ", as.character(crs(biovar, proj = T, describe = T)), "\n",
            "Region", as.character(crs(region, proj = T, describe = T)), "\n",
            "Making them identical...", "\n"
        )
      } else crs(region) = crs(biovars)
    }
    
    #Crop WorldClim data to region
    cat("Cropping and masking ", names(biovar), "\n")
    crop = crop(biovar, region)
    masked = mask(crop, region)
    biovarsMask[[length(biovarsMask) + 1]] = masked
  }
  
  if (plot == T) {
    for (i in 1:length(biovarsMask)) {
      if(verbose == T) cat("Plotting Bio variable ",as.character(i), "\n")
      # Plot each raster
      plot(biovarsMask[[i]], main = paste("Biovar", i))
    }
  }
  else cat("Plotting skipped \n")
  
  cat(cc$aquamarine("Cropped all biovars to region \n"))
  return(biovarsMask)
}
