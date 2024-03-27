wc_to_region <- function(biovars, region, projection, show_plot = F, verbose = T) {
  
  vebcat("Initiating WorldClim to region crop protocol.", color = "funInit")

  biovarsMask <- list()
  
  for (i in 1:terra::nlyr(biovars)) {
    biovar <- biovars[[i]]
    if (!isTRUE(identical(crs(biovar), crs(region)))) {
      vebcat(colcat("Crs not identical: \n", color = "nonFatalError"), 
             "Biovars: ", as.character(crs(biovar, proj = T, describe = T)), 
             "\n", "Region", as.character(crs(region, proj = T, describe = T)), 
             "\n", "Making them identical...", veb = verbose)
        
        if (projection == "longlat") {
          
          vebcat("Projecting to longlat.", veb = verbose)
          region <- project(region, crs(biovar))
          
        } else if (projection == "laea") {
          
          vebcat("Projecting to laea for Layer", highcat(i), "/", highcat(length(1:terra::nlyr(biovars))), veb = verbose)
          
          biovar <- terra::project(biovar, laea_crs)
          
        } else {
          stop("Missing or wrong use of projection parameter.")
        }
      
      vebcat("Crs made identical:", color = "proSuccess", veb = verbose)
      vebcat(
        "Biovars: ", as.character(crs(biovar, proj = T)), 
        "\n","Region", as.character(crs(region, proj = T)),
        veb = verbose
      )
      
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
    vebcat("Cropping and masking", highcat(sub("wc2.1_2.5m_", "", as.character(names(biovar)))), veb = verbose)
    crop <- crop(biovar, region)
    masked <- mask(crop, region)
    biovarsMask <- c(biovarsMask, list(masked))
  }
  
  biovarsMask <- terra::rast(biovarsMask)

  if (show_plot == T) {
    for (i in 1:terra::nlyr(biovarsMask)) {
      vebcat(paste0("Plotting bio_", as.character(i)), veb = verbose)
      # Plot each raster
      plot(biovarsMask[[i]], main = paste("Biovar", i))
    }
  } else {
    vebcat("Plotting skipped.", veb = verbose)
  }

  vebcat("WorldClim to region cropping protocol completed successfully", color = "funSuccess")

  
  names(biovars) <- sub("wc2.1_2.5m_", "", (names(biovars)))
  varnames(biovars) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))
  
  return(biovarsMask)

}
