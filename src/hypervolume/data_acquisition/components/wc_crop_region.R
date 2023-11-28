wc_to_region <- function(region = NULL, var = "bio", res = 2.5, path = "./resources/region", version = "2.1", projection, plot_it = F, verbose = T) {
  
  biovars <- get_wc_data()
  
  if (is.null(region)) {
    
    cat("Getting world data, skipping crop phase. \n")
    
    if (plot_it == T) {
      for (i in 1:length(biovars)) {
        cat(paste0("Plotting bio_", as.character(i)), "\n")
        # Plot each raster
        plot(biovars[[i]], main = paste("Biovar", i))
      }
    } else {
      cat("Plotting skipped \n")
    }

    biovars <- terra::rast(biovars)
    
    
    names(biovars) <- sub("wc2.1_2.5m_", "", (names(biovars)))
    varnames(biovars) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))

    return(biovars)
  } else {
    
    cat("Starting cropping sequence \n")

    biovarsMask <- list()

    for (biovar in biovars) {
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
            cat("Projecting to polar. \n")
            biovar <- terra::project(biovar, crs(region))
          } else {
            stop("Missing or wrong use of projection parameter.")
          }

          cat(
            green("Crs fixed: \n"),
            "Biovars: ", as.character(crs(biovar, proj = T, describe = T)), "\n",
            "Region", as.character(crs(region, proj = T, describe = T)), "\n",
            "Making them identical...", "\n"
          )
        } else {
          if (projection == "longlat") {
            biovar <- terra::project(biovar, crs(region))
          } else if (projection == "polar") {
            crs(region) <- crs(biovar)
          } else {
            stop("Missing or wrong use of projection parameter.")
          }
        }
      }

      # Crop WorldClim data to region
      if (verbose == T) cat("Cropping and masking", sub("wc2.1_2.5m_", "", as.character(names(biovar))), "\n")
      crop <- crop(biovar, region)
      masked <- mask(crop, region)
      biovarsMask[[length(biovarsMask) + 1]] <- masked
    }

    if (plot_it == T) {
      for (i in 1:length(biovarsMask)) {
        cat(paste0("Plotting bio_", as.character(i)), "\n")
        # Plot each raster
        plot(biovarsMask[[i]], main = paste("Biovar", i))
      }
    } else {
      cat("Plotting skipped \n")
    }

    cat(cc$aquamarine("Cropped all biovars to region \n"))

    biovars <- terra::rast(biovarsMask)

    
    names(biovars) <- sub("wc2.1_2.5m_", "", (names(biovars)))
    varnames(biovars) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))
    
    return(biovars)
  }
}
