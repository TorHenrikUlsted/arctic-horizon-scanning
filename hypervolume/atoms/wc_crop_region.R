wc_to_region = function(region, var = "bio", res = 2.5, path = "./resources", version = "2.1", plot = T, verbose= T) {
  biovars = list()
  # Check if all 19 biovariables are present, if not install
  for (i in 1:19) {
    file_name = paste0(path,"/wc2.1_",res,"m/","wc2.1_",res,"m_",var,"_",i,".tif")
    
    if(!file.exists(file_name)) {
      cat("WorldClim data does not exist, downloading... \n")
      r = worldclim_global(var, res, path, version)
    } else r = rast(file_name)
    
    biovars[[length(biovars) + 1]] = r
  }
  
  #print(biovars)
  
  biovarsMask = list()

  for (biovar in biovars) {
    #Crop WorldClim data to Arctic CAVM
    cat("Cropping and masking ", names(biovar), "\n")
    crop = crop(biovar, region)
    
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
