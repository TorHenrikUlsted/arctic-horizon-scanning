wc_to_region <- function(biovars, region, projection, show_plot = F, verbose = T) {
  
  cat("Initiating crop sequence \n")

  biovarsMask <- list()
  
  for (i in 1:terra::nlyr(biovars)) {
    biovar <- biovars[[i]]
    if (!isTRUE(identical(crs(biovar), crs(region)))) {
      if (verbose == T) {
        cat(
          red("Crs not identical: \n"),
          "Biovars: ", as.character(crs(biovar, proj = T, describe = T)), "\n",
          "Region", as.character(crs(region, proj = T, describe = T)), "\n",
          "Making them identical...", "\n"
        )
        
        if (projection == "longlat") {
          
          cat("Projecting to longlat. \n")
          crs(region) <- crs(biovar)
          
        } else if (projection == "laea") {
          
          cat("Projecting to laea for Layer", cc$lightSteelBlue(i), "/", cc$lightSteelBlue(length(1:terra::nlyr(biovars))),  "\n")
          
          biovar <- terra::project(biovar, cavm_crs)
          
        } else {
          stop("Missing or wrong use of projection parameter.")
        }

        cat(
          green("Crs made identical: \n"),
          "Biovars: ", as.character(crs(biovar, proj = T, describe = T)), "\n",
          "Region", as.character(crs(region, proj = T, describe = T)), "\n"
        )
      } else {
        if (projection == "longlat") {
          crs(region) <- crs(biovar)
        } else if (projection == "laea") {
          biovar <- terra::project(biovar, cavm_crs)
        } else {
          stop("Missing or wrong use of projection parameter.")
        }
      }
    }

    # Crop WorldClim data to region
    if (verbose == T) cat("Cropping and masking", cc$lightSteelBlue(sub("wc2.1_2.5m_", "", as.character(names(biovar)))), "\n")
    crop <- crop(biovar, region)
    masked <- mask(crop, region)
    biovarsMask <- c(biovarsMask, list(masked))
  }
  
  biovarsMask <- terra::rast(biovarsMask)

  if (show_plot == T) {
    for (i in 1:terra::nlyr(biovarsMask)) {
      cat(paste0("Plotting bio_", as.character(i)), "\n")
      # Plot each raster
      plot(biovarsMask[[i]], main = paste("Biovar", i))
    }
  } else {
    cat("Plotting skipped \n")
  }

  cat(cc$aquamarine("Cropped all biovars to region \n"))

  
  names(biovars) <- sub("wc2.1_2.5m_", "", (names(biovars)))
  varnames(biovars) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))
  
  return(biovarsMask)

}
