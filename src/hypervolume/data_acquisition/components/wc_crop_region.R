wc_to_region <- function(biovars, region, projection, show_plot = F, verbose = T) {
  
  if (verbose) cat(blue("Initiating WorldClim to region crop protocol. \n"))

  biovarsMask <- list()
  
  for (i in 1:terra::nlyr(biovars)) {
    biovar <- biovars[[i]]
    if (!isTRUE(identical(crs(biovar), crs(region)))) {
      if (verbose) cat(red("Crs not identical: \n"), "Biovars: ", as.character(crs(biovar, proj = T, describe = T)), "\n", "Region", as.character(crs(region, proj = T, describe = T)), "\n", "Making them identical...", "\n")
        
        if (projection == "longlat") {
          
          if (verbose) cat("Projecting to longlat. \n")
          region <- project(region, crs(biovar))
          
        } else if (projection == "laea") {
          
          if (verbose) cat("Projecting to laea for Layer", cc$lightSteelBlue(i), "/", cc$lightSteelBlue(length(1:terra::nlyr(biovars))),  "\n")
          
          biovar <- terra::project(biovar, laea_crs)
          
        } else {
          stop("Missing or wrong use of projection parameter.")
        }

      if (verbose) cat(green("Crs made identical: \n"),"Biovars: ", as.character(crs(biovar, proj = T, describe = T)), "\n","Region", as.character(crs(region, proj = T, describe = T)), "\n")
      } else {
        if (projection == "longlat") {
          region <- project(region, crs(biovar))
        } else if (projection == "laea") {
          biovar <- terra::project(biovar, laea_crs)
        } else {
          stop("Missing or wrong use of projection parameter.")
        }
      }

    # Crop WorldClim data to region
    if (verbose) cat("Cropping and masking", cc$lightSteelBlue(sub("wc2.1_2.5m_", "", as.character(names(biovar)))), "\n")
    crop <- crop(biovar, region)
    masked <- mask(crop, region)
    biovarsMask <- c(biovarsMask, list(masked))
  }
  
  biovarsMask <- terra::rast(biovarsMask)

  if (show_plot == T) {
    for (i in 1:terra::nlyr(biovarsMask)) {
      if (verbose) cat(paste0("Plotting bio_", as.character(i)), "\n")
      # Plot each raster
      plot(biovarsMask[[i]], main = paste("Biovar", i))
    }
  } else {
    if (verbose) cat("Plotting skipped \n")
  }

  if (verbose) cat(cc$aquamarine("WorldClim to region cropping protocol completed successfully \n"))

  
  names(biovars) <- sub("wc2.1_2.5m_", "", (names(biovars)))
  varnames(biovars) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))
  
  return(biovarsMask)

}
