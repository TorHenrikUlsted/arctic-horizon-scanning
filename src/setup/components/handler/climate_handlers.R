get_wc_data <- function(var = "bio", res = 2.5, path = "./resources/region", version = "2.1", show.plot = F, verbose = F) {
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

scale_biovars <- function(biovars, verbose = F) {
  name <- deparse(substitute(biovars))
  vebcat("Scaling biovariables", name, color = "funInit")

  save_dir <- "./outputs/setup/region"
  create_dir_if(save_dir)

  save_path <- paste0(save_dir, "/", name, "_scaled.tif")
  vebcat("save_path:", save_path, veb = verbose)

  # Check if the scaled data already exists
  if (file.exists(save_path)) {
    catn("Scaled data found, loading from file.")
    scaled_biovars <- rast(save_path)
    vebcat("Scaling biovariables completed successfully.", color = "funSuccess")
    return(scaled_biovars)
  }

  scaled_biovars <- biovars
  # Scale the dimensions
  for (i in 1:nlyr(scaled_biovars)) {
    catn("Scaling", highcat(names(scaled_biovars)[i]))
    # Extract the i-th layer
    layer <- biovars[[i]]

    # Scale the layer
    layer_scaled <- app(layer, fun = scale)

    # Replace the i-th layer of r_scaled with the scaled layer
    scaled_biovars[[i]] <- layer_scaled
  }

  vebcat("Renaming layers.", veb = verbose)
  names(scaled_biovars) <- names(biovars)

  # Save the scaled data
  writeRaster(scaled_biovars, save_path)
  catn(name, "scaled data saved to file:", colcat(save_path, color = "output"))

  vebcat("Scaling biovariables completed successfully.", color = "funSuccess")
  return(scaled_biovars)
}

wc_to_region <- function(biovars, shapefile, projection, show.plot = FALSE, verbose = FALSE) {
  vebcat("Initiating WorldClim to region crop protocol.", color = "funInit")

  save_path <- paste0("./outputs/setup/region/biovars-region-", projection, ".tif")

  if (file.exists(save_path)) {
    catn("Biovars already cropped to region, reading file..")
    biovarsMask <- rast(save_path)
  } else {
    region <- load_region(shapefile, verbose = verbose)

    biovarsMask <- list()

    region_crop_timer <- start_timer("crop_region")

    for (i in 1:terra::nlyr(biovars)) {
      biovar <- biovars[[i]]
      if (!isTRUE(identical(crs(biovar), crs(region)))) {
        vebcat(colcat("Crs not identical: \n", color = "nonFatalError"),
          "Biovars: ", as.character(crs(biovar, proj = T, describe = T)),
          "\n", "Region", as.character(crs(region, proj = T, describe = T)),
          "\n", "Making them identical...",
          veb = verbose
        )

        if (projection == "longlat") {
          vebcat("Projecting to longlat.", veb = verbose)
          region <- project(region, crs(biovar))
        } else if (projection == "laea") {
          vebcat("Projecting to laea for Layer", highcat(i), "/", highcat(length(1:terra::nlyr(biovars))), veb = verbose)

          biovar <- terra::project(biovar, config$projection$laea)
        } else {
          stop("Missing or wrong use of projection parameter.")
        }

        vebcat("Crs made identical:", color = "proSuccess", veb = verbose)
        vebcat(
          "Biovars: ", as.character(crs(biovar, proj = T)),
          "\n", "Region", as.character(crs(region, proj = T)),
          veb = verbose
        )
      } else {
        if (projection == "longlat") {
          region <- project(region, crs(biovar))
        } else if (projection == "laea") {
          biovar <- terra::project(biovar, config$projection$laea)
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

    if (show.plot) {
      for (i in 1:terra::nlyr(biovarsMask)) {
        vebcat(paste0("Plotting bio_", as.character(i)), veb = verbose)
        # Plot each raster
        plot(biovarsMask[[i]], main = paste("Biovar", i))
      }
    } else {
      vebcat("Plotting skipped.", veb = verbose)
    }

    names(biovarsMask) <- sub("wc2.1_2.5m_", "", (names(biovars)))
    varnames(biovarsMask) <- sub("wc2.1_2.5m_", "", (varnames(biovars)))

    catn("Writing out region raster to:", colcat(save_path, color = "output"))

    writeRaster(biovarsMask, save_path)

    end_timer(region_crop_timer)
  }

  vebcat("WorldClim to region cropping protocol completed successfully", color = "funSuccess")

  return(biovarsMask)
}
