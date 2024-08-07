get_climate_database <- function() {
  worldclim_args <- c("wc", "worldclim", "world-clim", "world_clim")
  chelsa_args <- c("c", "ch", "chelsa")
  
  if (tolower(config$run$climate$database) %in% tolower(worldclim_args)) {
    database <- "worldclim"
  } else if (tolower(config$run$climate$database) %in% tolower(chelsa_args)) {
      database <- "chelsa"
  } else {
    vebcat("Database not found. Edit run file to use 'worldclim' or 'chelsa' database.", color = "fatalError")
    stop("Edit run parameter 'climate.database'")
  }
  
  return(database)
}

download_climate_data <- function(database = "worldclim", verbose = FALSE) {
  if (database == "worldclim") {
    invisible(worldclim_global(
      var = config$climate$database$worldclim$var, 
      res = config$climate$database$worldclim$res,
      path = "./resources", 
      version = config$climate$database$worldclim$version
    ))
  } else if (database == "chelsa") {
    chelsa_download(
      dataset = config$climate$database$chelsa$dataset, 
      variable = config$climate$database$chelsa$var, 
      period = config$climate$database$chelsa$period, 
      model = config$climate$database$chelsa$model, 
      ssp = config$climate$database$chelsa$ssp, 
      dir.out = "./resources/climate/chelsa", 
      version = config$climate$database$chelsa$version, 
      verbose = verbose
    )
  } else {
    vebcat("Error: database not found.", color = "fatalError")
    stop("Edit config in './src/utils/config.R'")
  }
}

handle_climate_filename <- function(database = "worldclim", verbose = FALSE) {
  path = "./resources"
  if (database == "worldclim") {
    
    dir_path <- paste0(path, "/climate/wc", config$climate$database$worldclim$version, "_", config$climate$database$worldclim$res, "m")
    
    file_names <- names(worldclim_global(config$climate$database$worldclim$var, config$climate$database$worldclim$res, "./resources"))
    
    file_names <- paste0(dir_path, "/", file_names, ".tif")
    
  } else if (database == "chelsa") {
    dir_path <- paste0(path, "/climate/chelsa/", var)
    file_names <- handle_chelsa_filenames(
      dataset = config$climate$database$chelsa$dataset, 
      variable = config$climate$database$chelsa$var, 
      model = config$climate$database$chelsa$model, 
      ssp = config$climate$database$chelsa$ssp, 
      period = config$climate$database$chelsa$period, 
      version = config$climate$database$chelsa$version, 
      verbose = verbose
    )
    
    file_names <- paste0(dir_path, "/", file_names)
    
  } else {
    vebcat("Cannot find database, check config settings.", color = "fatalError")
    stop("Choose another database.")
  }
  
  return(file_names)
} 

load_climate_data <- function(database, show.plot = FALSE, verbose = FALSE) {
  biovars <- list()
  
  database <- get_climate_database()
  
  download_climate_data(database, verbose = FALSE)
  
  file_names <- handle_climate_filename(database, verbose = FALSE)

  biovars <- terra::rast(file_names)
  
  names(biovars) <- paste0("bio", seq(file_names))
  varnames(biovars) <- paste0("bio", seq(file_names))
  
  return(biovars)
}

handle_biovar_saves <- function() {
  database <- get_climate_database()
  
  save_dir <- paste0("./outputs/setup/region/", database)
  
  if (database == "worldclim") {
    save_dir <- paste0(save_dir, "_", config$climate$database$worldclim$version, "/", config$climate$database$worldclim$res, "/", config$climate$database$worldclim$period)
  } else if (database == "chelsa") {
    if (!is.null(config$climate$database$chelsa$model) && !is.null(config$climate$database$chelsa$ssp)) {
      save_dir <- paste0(save_dir, "_0.5/", config$climate$database$chelsa$period, "/", config$climate$database$chelsa$model, "_", config$climate$database$chelsa$ssp)
    }
    save_dir <- paste0(save_dir, "_0.5/", config$climate$database$chelsa$period)
  }
  
  return(save_dir)
}

scale_biovars <- function(biovars, verbose = FALSE) {
  name <- deparse(substitute(biovars))
  vebcat("Scaling biovariables", name, color = "funInit")
  
  scale_timer <- start_timer("climate-scaling-timer")
  
  save_dir <- handle_biovar_saves()
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
  
  vebcat("Creating scale template")
  # Create a template for the output raster
  
  print(biovars)
  print(nlyr(biovars))
  
  # Start the writing process
  filecon  <- writeStart(biovars, filename = save_path, overwrite = FALSE)
  
  print(filecon)
  
  # Process and write each layer
  for (i in 1:nlyr(biovars)) {
    catn("Scaling", highcat(names(biovars)[i]))
    
    # Extract the entire layer
    layer <- biovars[[i]]
    print(layer)
    
    # Scale the layer
    layer_scaled <- scale(layer)
    print(layer_scaled)
    
    readValues(layer_scaled)
    
    # Write the scaled layer
    #writeValues(layer_scaled, filecon, start = layer, nrows = filecon$nrows)
  }
  
  # Finish the writing process
  writeStop(filecon)
  
  catn(name, "scaled data saved to file:", colcat(save_path, color = "output"))
  
  scaled_biovars <- rast(save_path)

  end_timer(scale_timer)

  vebcat("Scaling biovariables completed successfully.", color = "funSuccess")
  print(scaled_biovars)
  return(scaled_biovars)
}

climate_to_region <- function(biovars, shapefile, projection, show.plot = FALSE, verbose = FALSE) {
  vebcat("Initiating WorldClim to region crop protocol.", color = "funInit")

  save_dir <- handle_biovar_saves()
  create_dir_if(save_dir)
  
  save_path <- paste0(save_dir, "/biovars-region-", projection, ".tif")

  if (file.exists(save_path)) {
    catn("Biovars already cropped to region, reading file..")
    biovarsMask <- rast(save_path)
  } else {
    region <- load_region(shapefile, verbose = verbose)

    region_crop_timer <- start_timer("crop_region")
    
    biovar_template <- crop(biovars[[1]], region)
    biovar_template <- mask(biovar_template, region)
    
    out <- writeStart(biovar_template, save_path, overwrite=TRUE, n=nlyr(biovars))
    
    writeValues(out, values(biovar_template), start=1, nrows=nrow(biovar_template))

    for (i in 2:terra::nlyr(biovars)) {
      biovar <- biovars[[i]]
      
      check_crs(biovar, config$projection$crs$longlat, "bilinear", verbose)

      # Crop climate data to region
      vebcat("Cropping and masking", highcat(names(biovar)), veb = verbose)
      crop <- crop(biovar, region)
      masked <- mask(crop, region)
      
      names(biovarsMask) <- paste0("bio", i)
      varnames(biovarsMask) <- paste0("bio", i)
      
      writeValues(out, values(masked), start=1, nrows=nrow(masked))
    }
      
    out <- writeStop(out)

    biovarsMask <- terra::rast(save_path)

    if (show.plot) {
      for (i in 1:terra::nlyr(biovarsMask)) {
        vebcat(paste0("Plotting bio_", as.character(i)), veb = verbose)
        # Plot each raster
        plot(biovarsMask[[i]], main = paste("Biovar", i))
      }
    } else {
      vebcat("Plotting skipped.", veb = verbose)
    }

      end_timer(region_crop_timer)
    }
    
  catn("Raster written to:", colcat(save_path, color = "output"))

  vebcat("WorldClim to region cropping protocol completed successfully", color = "funSuccess")

  return(biovarsMask)
}
