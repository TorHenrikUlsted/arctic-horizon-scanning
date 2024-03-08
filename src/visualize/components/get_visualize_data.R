get_stats_data <- function(cleaned_data, hv.dir, hv.method, verbose = F, warn, err) {
  cat(blue("Initializing stats data. \n"))
  # Define directories
  vis_stats_dir <- "./outputs/visualize/stats"
  
  if (!file.exists("./outputs/visualize/stats/included-species.csv") || !file.exists("./outputs/visualize/stats/excluded-species.csv") ) {
    stats_src <- paste0(hv.dir, "/stats/box-stats.csv")
    sp_stats <- cleaned_data
    
    # Get clean species csv
    excluded_sp <- sp_stats %>% filter(excluded == TRUE)
    
    fwrite(excluded_sp, paste0(vis_stats_dir, "/uniq-excluded-species.csv"), bom = TRUE)
    
    included_sp <- sp_stats %>% filter(excluded == FALSE)
    
    fwrite(included_sp, paste0(vis_stats_dir, "/uniq-included-species.csv"), bom = TRUE)
    
    # Get the region for the species
    glonaf_wfo_one <- fread("./resources/synonym-checked/glonaf-species-wfo-one.csv")
    
    gwo <- glonaf_wfo_one %>%
      select(rawName.ORIG) %>%
      # Get the refined scientificNames as used in the hypervolume method
      mutate(species = apply(glonaf_wfo_one[, c("genus", "specificEpithet", "infraspecificEpithet"), drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " ")))
    
    gwo$species <- trimws(gwo$species)
    # Swap space with line to match the sp-stats
    gwo$species <- gsub(" ", "-", gwo$species)
    
    colnames(gwo) <- c("origName", "species")
    
    # Remove duplicates
    gwo <- gwo[!duplicated(gwo$species)]
    
    # Match the names with the sp_stats
    gwo_stats <- sp_stats %>%
      left_join(gwo, by = "species")
    
    # Get the region id
    glonaf_orig_df <- fread("./resources/data-raw/glonaf.csv")
    
    godf <- glonaf_orig_df %>%
      select(standardized_name, region_id)
    
    colnames(godf) <- c("origName", "regionId")
    
    glonaf_orig_region <- fread("./resources/region/glonaf-region-desc.csv")
    
    gor <- glonaf_orig_region %>% select(region_id, country_ISO, country)
    
    colnames(gor) <- c("regionId", "countryIso", "country")
    
    # Match the gwo_stats and godf --- Some will have multiple region ids
    matched_glonaf <- merge(gwo_stats, godf, by = "origName")
    
    if (verbose) cat("Number of rows when appending region ids:", cc$lightSteelBlue(nrow(matched_glonaf)), "\n")
    
    # Remove duplicates based on combined region id and standardized name
    matched_dups <- nrow(matched_glonaf[duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ])
    if (verbose) cat("Number of duplicates after appending region id:", cc$lightSteelBlue(matched_dups), "\n")
    
    if (matched_dups > 0) {
      matched_glonaf <- matched_glonaf[!duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ]
      
      if (verbose) cat("Number of duplicates after removal:", cc$lightSteelBlue(nrow(matched_glonaf[duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ])))
      
      if (verbose) cat("Number of rows after finishing appending region ids:", cc$lightSteelBlue(nrow(matched_glonaf)), "\n")
    }
    
    # Append region info
    matched_stats <- merge(matched_glonaf, gor, by = "regionId")
    
    matched_stats <- matched_stats %>%
      select(species, origName, observations, dimensions, samplesPerPoint, randomPoints, excluded, jaccard, sorensen, fracVolumeSpecies, fracVolumeRegion, overlapRegion, includedOverlap, regionId, countryIso, country)
    
    # Seperate into included and excluded species
    
    included_sp <- matched_stats %>% filter(excluded == FALSE)
    
    if (verbose) cat("Number of included species:", cc$lightSteelBlue(nrow(included_sp)), "\n")
    
    create_dir_if("./outputs/visualize/stats")
    
    fwrite(included_sp, paste0(vis_stats_dir, "/included-species.csv"), bom = TRUE)
    
    excluded_sp <- matched_stats %>% filter(excluded == TRUE)
    
    if (verbose) cat("Number of included species:", cc$lightSteelBlue(nrow(excluded_sp)), "\n")
    
    fwrite(excluded_sp, paste0(vis_stats_dir, "/excluded-species.csv"), bom = TRUE)
    
  } else {
    included_sp <- fread(paste0(vis_stats_dir, "/included-species.csv"))
    
    excluded_sp <- fread(paste0(vis_stats_dir, "/excluded-species.csv"))
  }
  
  cat(cc$lightGreen("Stats successfully initialised.\n"))
  
  return(included_sp)
}

get_inclusion_data <- function(region, hv.dir, hv.method, hv.project.method, verbose = F, warn, err) {
  cat(blue("Initializing visulasation data. \n"))
  
  save_file <- "./outputs/visualize/logs/get-visualize-data/cell-regions-dt.csv"
  
  if (file.exists("./outputs/visualize/logs/get-visualize-data/cell-regions-dt.csv")) {
    cell_regions_dt <- fread(save_file)
  } else {
    sp_dirs <- list.dirs(paste0(hv.dir, "/projections/", hv.method))
    # Remove the directory name
    sp_dirs <- sp_dirs[-1]
    
    if (verbose) {
      cat("Sample directory names: \n")
      print(head(sp_dirs, 3))
    }
    
    # Initialize  data table
    cat("Initializing data table.\n")
    sp_rast <- terra::rast(paste0(sp_dirs[[1]], "/", hv.project.method, ".tif"))
    names(sp_rast) <- basename(sp_dirs[[1]])
    
    if (!identical(ext(sp_rast), ext(region))) {
      cat("Cropping region extent to match species.\n")
      region <- crop(region, ext(sp_rast))
    }
    
    region_rast_dt <- terra::extract(sp_rast, region, cells = TRUE, na.rm = FALSE)
    region_rast_dt <- data.table(region_rast_dt)
    
    region_dt <- as.data.table(region)
    region_dt[, ID := .I]
    
    cell_regions_dt <- merge(region_rast_dt, region_dt[, .(ID, FLOREG, country, floristicProvince)], by = "ID")
    
    cell_regions_dt <- cell_regions_dt[, -2]
    # Merge cell regions with the raster cells
    sp_dt <- terra::extract(sp_rast, ext(sp_rast), cells = TRUE)
    sp_dt <- data.table(sp_dt)
    names(sp_dt) <- c("cell", "richness")
    
    cell_regions_dt <- merge(sp_dt, cell_regions_dt[, .(cell, FLOREG, country, floristicProvince)], by = "cell", all = TRUE)
    
    cat("Table sample:\n")
    print(head(cell_regions_dt, 2))
    
    #
    for (i in 2:length(sp_dirs)) {
      species <- sp_dirs[[i]]
      
      cat("\rCalculating richness in raster:", i, "/", length(sp_dirs))
      
      flush.console()
      
      sp_file <- paste0(species, "/", hv.project.method, ".tif")
      sp_name <- basename(species)
      
      sp_rast <- terra::rast(sp_file)
      names(sp_rast) <- sp_name
      #t <- global(sp_rast, fun = sum, na.rm=TRUE) # Sum all values in the raster
      
      sp_rast_dt <- terra::extract(sp_rast, ext(sp_rast), cells = TRUE)
      sp_rast_dt <- data.table(sp_rast_dt)
      names(sp_rast_dt) <- c("cell", "value")
      
      # Merge to sum with the richness column
      cell_regions_dt <- merge(cell_regions_dt, sp_rast_dt, by = "cell", all = TRUE)
      
      # Sum richness with new value, only if both are not NA, else keep richness value
      cell_regions_dt[, richness := ifelse(!is.na(richness) & !is.na(value), richness + value, richness), by = cell]
      cell_regions_dt[, value := NULL]
    }; cat("\n")
    
    create_dir_if("./outputs/visualize/logs/get-visualize-data")
    
    fwrite(cell_regions_dt, save_file, bom = TRUE)
  }

  cat(cc$lightGreen("visualization data successfully initialised.\n"))

  return(cell_regions_dt)
}

get_projection_max <- function(hv.dir, hv.method, hv.project.method  = "probability", n = NULL) {
  
  cat(blue("Acquiring max values for all probability rasters. \n"))
  
  log_file <- paste0("./outputs/visualize/logs/", hv.project.method, "-max-values.csv")
  
  if (file.exists(log_file)) {
    cat("Found max value table:", log_file, "\n")
    
    max_values <- fread(log_file)
    
  } else {
    cat("No previous table found, acquiring max values for", hv.project.method, "method.\n")
    
    sp_dirs <- list.dirs(paste0(hv.dir, "/projections/", hv.method))
    
    # Remove the directory name
    sp_dirs <- sp_dirs[-1]
    
    max_values <- data.table()
    
    for (i in seq_along(sp_dirs)) {
      sp_dir <- sp_dirs[[i]]
      sp_name <- basename(sp_dir)
      sp_filename <- paste0(sp_dir, "/", hv.project.method, ".tif")
      
      cat("\rAcquiring max values", i, "/", length(sp_dirs))
       
      sp_rast <- terra::rast(sp_filename)
      sp_max <- where.max(sp_rast)
      sp_max <- data.table(sp_max)
      sp_max <- sp_max[, -(1:2)]
      
      sp_dt <- sp_max[, .(maxValue = value, species = sp_name, filename = sp_filename)]
      
      max_values <- rbind(max_values, sp_dt)
      
      max_values <- max_values[maxValue != 0, ]
    }; cat("\n")
    
    max_values <- max_values[maxValue != 0, ]
    
    max_values <- max_values[order(-max_values$maxValue), ]
    
    fwrite(max_values, log_file, bom = TRUE)
  }
  
  if (is.null(n)) {
    cat("n is not used, returning the whole table.\n")
  } else {
    max_values <- max_values[1:n, ]
  }
  
  cat(cc$lightGreen("Successfully acquired max values for all probability rasters.\n"))
  
  return(max_values)
}

get_prob_max <- function(spec.filename) {
  sp_name <- basename(dirname(spec.filename))
  
  sp_rast <- terra::rast(spec.filename)
  sp_max <- where.max(sp_rast)
  sp_max <- sp_max[sp_max != 0]
  sp_max <- data.table(sp_max)
  sp_max <- sp_max[, -(1:2)]
  
  sp_dt <- sp_max[, .(maxValue = value, species = sp_name, filename = spec.filename)]
  
  return(sp_dt)
}

get_inclusion_coverage <- function(spec.filename) {
  sp_name <- basename(dirname(spec.filename))
  
  sp_rast <- terra::rast(spec.filename)
  
  sp_coverage <- terra::global(sp_rast, fun = sum, na.rm = TRUE)
  
  sp_dt <- data.table(coverage = sp_coverage, species = sp_name, filename = spec.filename)
  
  return(sp_dt)
}

stack_projections <- function(filenames, projection, projection.method, binary = FALSE, verbose = FALSE) {
  
  cat(blue("Stacking probability rasters. \n"))
  # Initialize raster stack with the first layer
  sp_name <- filenames[['species']][[1]]
  sp_filename <- filenames[['filename']][[1]]
  sp_rast <- terra::rast(sp_filename)
  names(sp_rast) <- sp_name
  
  if(!identical(crs(sp_rast, proj = TRUE), crs(projection, proj = TRUE))) {
    cat("Reprojecting to", as.character(crs(projection, proj = TRUE)), "\n")
    stack <- terra::project(sp_rast, projection, method = projection.method)
  } else {
    stack <- sp_rast
  }
  
  for (i in 2:nrow(filenames)) {
    sp_name <- filenames[['species']][[i]]
    sp_filename <- filenames[['filename']][[i]]
    
    cat("\rStacking raster", i, "/", nrow(filenames))
    
    sp_rast <- terra::rast(sp_filename)
    names(sp_rast) <- sp_name
    
    if(!identical(crs(sp_rast, proj = TRUE), crs(projection, proj = TRUE))) {
      if (verbose) cat("Reprojecting to", as.character(crs(projection, proj = TRUE)), "\n")
      sp_rast <- terra::project(sp_rast, projection, method = projection.method)
    }
    
    stack <- c(stack, sp_rast)

  }; cat("\n")
  
  # Convert values to 0 and 1 if binary
  if (binary) stack <- ceiling(stack)
  
  cat(cc$lightGreen("Successfully stacked probability rasters.\n"))
  
  return(stack)
}

convert_template_raster <- function(input.values, hv.dir, hv.method, hv.project.method, projection, projection.method = "near", verbose = F) {
  cat(blue("Acquiring template raster\n"))
  sp_dirs <- list.dirs(paste0(hv.dir, "/projections/", hv.method))
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]
  
  sp_filename <- paste0(sp_dirs[[1]], "/", hv.project.method, ".tif")
  
  template <- terra::rast(sp_filename)
  # Check length
  cat("Length of input / raster\n")
  cat(length(input.values), "/", ncell(template),"\n")
  
  terra::values(template) <- input.values
  
  template <- check_crs(template, projection = projection, projection.method = projection.method)
  
  return(template)
}

get_sp_names_region <- function(spec.filename) {
  
}

get_world_map <- function(projection, scale = "medium", pole = "north") {
  cat(blue("Setting up world Map.\n"))
  
  world_map <- ne_countries(scale = scale, returnclass = "sf")
  
  world_map <- vect(world_map)
  
  if (pole == "north") {
    equator_extent <- ext(-180, 180, 0, 90)
  } else if (pole == "south") {
    equator_extent <- ext(-180, 180, -90, 0)
  } else {
    stop("pole has to be either 'south' or 'north'.")
  }
  
  world_map <- crop(world_map, equator_extent)
  
  world_map <- project(world_map, projection)
  
  cat(cc$lightGreen("World Map setup completed successfully.\n"))
  
  return(world_map)
}