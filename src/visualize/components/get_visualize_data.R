get_stats_data <- function(vis.dir, hv.dir, hv.method, projection = "longlat", verbose = F, warn, err) {
  vebcat("Initializing stats data.", color = "funInit")
  # Define directories
  stats_dir <- paste0(vis.dir, "/stats")
  create_dir_if(stats_dir)
  
  stats_file <- paste0(hv.dir, "/", hv.method, "-sequence/stats/stats.csv")
  inc_sp_out <- paste0(stats_dir, "/included-species.csv")
  exc_sp_out <- paste0(stats_dir, "/excluded-species.csv")
  
  if (!file.exists(inc_sp_out) || !file.exists(exc_sp_out) ) {
    
    stats <- fread(stats_file)
    
    included_sp <- stats[excluded == FALSE, ]
    
    fwrite(included_sp, inc_sp_out, bom = TRUE)
    
    vebcat("Number of included species:", highcat(length(unique(included_sp$cleanName))), veb = verbose)
    
    excluded_sp <- stats[excluded == TRUE, ]
    
    fwrite(excluded_sp, exc_sp_out, bom = TRUE)
    
    vebcat("Number of excluded species:", highcat(length(unique(excluded_sp$cleanName))), veb = verbose)
    
    rm(stats)
    invisible(gc())
    
  } else {
    included_sp <- fread(inc_sp_out)
  }
  
  vebcat("Stats successfully initialised.", color = "funSuccess")
  
  return(included_sp)
}

get_inclusion_cell <- function(spec.filename, region = NULL, extra = NULL, verbose = FALSE) {
  
  init <- function(spec.filename, region, verbose) {
    vebcat("Initiating get inclusion cell.", color = "funInit")
    
    if (is.null(region)) {
      vebcat("Missing region.")
      return(NULL)
    }
    
    sp_rast <- load_sp_rast(spec.filename)
    
    if (!identical(ext(sp_rast), ext(region))) {
      catn("Cropping region extent to match species")
      region <- crop(region, ext(sp_rast))
    }
    
    catn("Extracting raster from region.")
    sp_cell_dt <- extract_raster_to_dt(sp_rast, region, value = "cellAbundance", cells = TRUE)
    vebprint(head(sp_cell_dt, 3), verbose)
    
    catn("Converting region to data table.")
    region_dt <- as.data.frame(region)
    region_dt <- as.data.table(region_dt)
    region_dt <- handle_region_dt(region_dt)
    region_dt[, ID := .I]
    
    vebprint(head(region_dt, 3), verbose)
    
    catn("Merging raster_dt and region_dt.")
    sp_regions_dt <- merge(sp_cell_dt, region_dt, by = "ID")
    sp_regions_dt[, ID := NULL]
    
    vebprint(head(sp_regions_dt, 3), verbose)
    
    # Merge cell regions with the raster cells
    catn("Extracting raster cells.")
    sp_dt <- extract_raster_to_dt(sp_rast, value = "toRemove", cells = TRUE)
    sp_dt[, toRemove := NULL]
    
    catn("Merging raster cells with region by cells.")
    cell_regions_dt <- merge(sp_dt, sp_regions_dt, by = "cell", all = TRUE)
    
    # st unique cell and if any have duplicates and are above 0 set to 1
    cell_regions_dt[, cellAbundance := ifelse(any(cellAbundance > 0), 1, cellAbundance), by = "cell"]
    
    cell_regions_dt <- cell_regions_dt[order(cell_regions_dt$cell)]
    
    vebcat("get inclusion cell initialized successfully.", color = "funSuccess")
    
    return(cell_regions_dt)
  }
  
  execute <- function(spec.filename, region, extra = NULL, verbose) {
    
    catn("Initiating data table.")
    sp_name <- basename(dirname(spec.filename[1]))
    vebcat(sp_name, veb = verbose)
    
    sp_rast <- load_sp_rast(spec.filename[1])
    
    res_dt <- extract_ext_to_dt(sp_rast, value = "valueSum")
    vebprint(res_dt, veb = verbose)
    
    for(i in 2:length(spec.filename)) {
      catn("Execution Iteration:", i)
      species <- spec.filename[i]
      
      sp_name <- basename(dirname(species))
      catn(sp_name)
      
      sp_rast <- load_sp_rast(species)
      
      vebcat("Extracting raster cells.", veb = verbose)
      sp_rast_dt <- extract_ext_to_dt(sp_rast, value = "value")
      
      vebprint(head(sp_rast_dt, 3), veb = verbose)
      
      vebcat("Merging with init results.", veb = verbose)
      # Merge to sum with the richness column
      cell_regions_dt <- merge(res_dt, sp_rast_dt, by = "cell", all = TRUE)
      
      vebprint(head(cell_regions_dt, 3), veb = verbose)
      
      vebcat("Calculating richness", veb = verbose)
      # Sum richness with new value, only if both are not NA, else keep richness value
      cell_regions_dt[, valueSum := ifelse(!is.na(valueSum) & !is.na(value), valueSum + value, ifelse(!is.na(valueSum), valueSum, value)), by = cell]
      
      cell_regions_dt[, value := NULL] # Remove the value column after summing to save ram and space
      res_dt <- cell_regions_dt
      
      vebprint(head(res_dt, 3), veb = verbose)
    }
    
    catn("Finished Executing parallel.")
    
    return(res_dt)
  }
  
  # Handle the parallel results will return a list of results
  process <- function(parallel.res, verbose) { # parallel.res is a list of data tables
    dt_region <- parallel.res[[1]]
    dt_list <- parallel.res[[2]]
    
    summed_dt <- dt_region
    
    vebprint(dt_region, veb = verbose)
    vebprint(dt_list, veb = verbose)
    
    for (i in 1:length(dt_list)) {
      cat("\rProcessing list", i, "/", length(dt_list))
      
      dt_cell <- dt_list[[i]]
      
      merged_dt <- merge(summed_dt, dt_cell, by = "cell")
      
      # Sum the richness
      merged_dt[, richness := ifelse(!is.na(richness) & !is.na(valueSum), richness + valueSum, ifelse(!is.na(richness), richness, valueSum)), by = cell]
      
      merged_dt[, valueSum := NULL]
      
      summed_dt <- merged_dt
    }; catn()
    
    
    return(summed_dt)
  }
  
  return(list(
    init = init,
    execute = execute,
    process = process
  ))
}

convert_template_raster <- function(input.values, hv.proj.dir, hv.method, hv.project.method, projection, projection.method = "near", out.dir, verbose = F) {
  
  out_file <- paste0(out.dir, "/", hv.project.method,"-hotspots.tif")
  
  if (file.exists(out_file)) {
    catn("Reading existing raster.")
    template <- terra::rast(out_file)
    return(template)
  } 
  
  vebcat("Converting raster", color = "funInit")
  
  sp_dirs <- list.dirs(hv.proj.dir)
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]
  
  sp_filename <- paste0(sp_dirs[[1]], "/", hv.project.method,".tif")
  
  template <- terra::rast(sp_filename)
  # Check length
  catn("Length of input / raster")
  catn(length(input.values), "/", ncell(template))
  
  terra::values(template) <- input.values
  
  template <- check_crs(template, projection = projection, projection.method = projection.method)
  
  catn("Writing raster to:", colcat(out_file, color = "output"))
  writeRaster(template, out_file)
  
  vebcat("Raster converted successfully", color = "funSuccess")
  
  return(template)
}

get_world_map <- function(projection, scale = "medium", pole = NULL) {
  vebcat("Setting up world Map.", color = "funInit")
  
  world_map <- ne_countries(scale = scale, returnclass = "sf")
  
  world_map <- vect(world_map)
  
  if (!is.null(pole)) {
    if (pole == "north") {
      equator_extent <- ext(-180, 180, 0, 90)
    } else if (pole == "south") {
      equator_extent <- ext(-180, 180, -90, 0)
    } else {
      stop("pole has to be either 'south' or 'north'.")
    }
    
    world_map <- crop(world_map, equator_extent)
  }
  
  world_map <- project(world_map, projection)
  
  vebcat("World Map setup completed successfully", color = "funSuccess")
  
  return(world_map)
}

get_inc_coverage <- function(spec.filename, region = NULL, extra = NULL, verbose = FALSE) {
  
  execute <- function(spec.filename, region, verbose) {
    sp_name <- basename(dirname(spec.filename))
    sp_name <- gsub("-", " ", sp_name)
    
    sp_rast <- terra::rast(spec.filename)
    
    catn("Extracting raster cells.")
    sp_dt <- terra::extract(sp_rast, ext(sp_rast), cells = TRUE)
    sp_dt <- as.data.frame(sp_dt)
    sp_dt <- as.data.table(sp_dt)
    names(sp_dt) <- c("cell", "value")
    
    # subset to only get values above 0
    sp_coverage <- sp_dt[value > 0]
    
    # Set all values to 1
    sp_coverage <- sp_coverage[, value := ifelse(value > 0, 1, value)]
    
    sp_coverage <- unique(sp_coverage, by = "cell")
    
    sp_coverage <- sum(sp_coverage$value)
    
    vebcat("Species included in", sp_coverage, "cells.")
    
    #sp_coverage <- terra::global(sp_rast, fun = sum, na.rm = TRUE)
    
    out_dt <- data.table(coverage = sp_coverage, species = sp_name, filename = spec.filename)
    
    print(out_dt)
    
    return(out_dt)
  }
  
  return(list(
    execute = execute
  ))
}

stack_projections <- function(filenames, projection, projection.method, out.dir, binary = FALSE, verbose = FALSE) {
  vebcat("Stacking rasters", color = "funInit")
  
  nodes_dir <- paste0(out.dir, "/nodes")
  rasters_dir <- paste0(out.dir, "/", gsub(".tif", "", basename(filenames)[1]) )
  
  create_dir_if(c(nodes_dir, rasters_dir))
  
  rast_out_names <- paste0(basename(dirname(filenames)), ".tif")
   
  raster_files <- list.files(rasters_dir, full.names = FALSE)
  
  skip <- FALSE
  
  if (all(rast_out_names %in% raster_files)) {
    catn("Found raster files, skipping to stacking.")
    skip <- TRUE
  }
  
 if (!skip) {
   max_cores <- calc_num_cores(
     ram.high = 5, 
     verbose =  FALSE
   )
   
   catn("Found", length(filenames), "files.")
   
   max_cores <- min(length(filenames), max_cores$high)
   
   catn("Creating cluster of", max_cores, "core(s).")
   
   cl <- makeCluster(max_cores)
   
   clusterEvalQ(cl, {
     source("./src/utils/utils.R")
   })
   
   filenames_dt <- data.table(
     order = seq_along(filenames),
     filename = filenames
   )
   
   # Get the exports that are not null
   export_vars <- c("filenames_dt", "projection", "projection.method", "out.dir", "binary", "verbose", "nodes_dir", "rasters_dir")
   
   vebcat("export_vars:", veb = verbose)
   vebprint(export_vars, veb = verbose)
   
   clusterExport(cl, export_vars, envir = environment())
   
   catn("Stacking rasters.")
   
   # Apply the function in parallel
   clusterApplyLB(cl, 1:length(filenames), function(i) {
     tryCatch({
       node_file <- paste0(nodes_dir, "/node", i)
       create_file_if(node_file)
       
       try(node_con <- file(node_file, "a"))
       sink(node_con, type = "output")
       sink(node_con, type = "message")
       
       file <- filenames_dt$filename[i]
       rast_out_names <- paste0(basename(dirname(file)), ".tif")
       name <- gsub("-", " ", basename(dirname(file)))
       
       out_file <- paste0(rasters_dir, "/", rast_out_names)
       
       catn("name:", name)
       catn("Input file:", file)
       catn("Output file:", out_file)
       catn("Input projection:", as.character(crs(projection, proj = TRUE)))
       catn("Projection method:", projection.method)
       catn()
       
       raster <- terra::rast(file)
       names(raster) <- name
       
       vebprint(raster, text = "Input raster:")
       
       if(!identical(crs(raster, proj = TRUE), crs(projection, proj = TRUE))) {
         catn()
         vebcat("Reprojecting to", as.character(crs(projection, proj = TRUE)))
         catn()
         try(raster <- terra::project(raster, projection, method = projection.method))
         vebprint(raster, text = "Output raster:")
       }
       
       # Convert values to 0 and 1 if binary
       if (binary) raster <- ceiling(raster)
       
       catn("Writing raster to:", colcat(out_file, color = "output"))
       terra::writeRaster(raster, out_file)
       
       sink(type = "message")
       sink(type = "output")
       close(node_con)
       
       invisible(gc())
     }, error = function(e) {
       try(node_con <- file(node_file, "a"))
       writeLines(paste0("Error when stacking rasters in parallel:\n", e), node_con)
       close(node_con)
       stopCluster(cl)
       closeAllConnections()
       stop()
     })
   })
   
   catn("finishing up.")
   
   stopCluster(cl)
   closeAllConnections()
 }
  
  catn("Listing files and stacking them.")
  raster_files <- list.files(rasters_dir, full.names = TRUE)
  
  # Order it by filenames order
  ordered <- raster_files[match(rast_out_names, basename(raster_files))]
  vebprint(ordered, verbose, "Ordered Rasters:")
  
  stack <- terra::rast(ordered[1])
  
  for (i in 2:length(ordered)) {
    cat("\rStacking", i, "/", length(ordered))
    
    file <- ordered[i]
    
    raster <- terra::rast(file)
    
    stack <- c(stack, raster)
  }; catn()
  
  vebcat("Successfully stacked probability rasters", color = "funSuccess")
  
  return(stack)
}

get_prob_stats <- function(spec.filename, region = NULL, extra, verbose = FALSE) {
  
  execute <- function(spec.filename, region, verbose) {
    sp_name <- basename(dirname(spec.filename))
    sp_name <- gsub("-", " ", sp_name)
    
    sp_rast <- terra::rast(spec.filename)
    
    sp_vals <- terra::values(sp_rast)
    
    #sp_vals <- sp_vals[!is.nan(sp_vals) & sp_vals != 0]
    
    sp_vals <- sp_vals[!is.nan(sp_vals)]
    
    if (length(sp_vals) == 0) {
      sp_min <- 0
      sp_mean <- 0
      sp_median <- 0
      sp_q1 <- 0
      sp_q2 <- 0
      sp_q3 <- 0
      sp_max <- 0
    } else {
      sp_min <- min(sp_vals, na.rm = TRUE)
      sp_mean <- mean(sp_vals, na.rm = TRUE)
      sp_median <- median(sp_vals, na.rm = TRUE)
      sp_q1 <- quantile(sp_vals, probs = 0.25, na.rm = TRUE)
      sp_q2 <- quantile(sp_vals, probs = 0.5, na.rm = TRUE)
      sp_q3 <- quantile(sp_vals, probs = 0.75, na.rm = TRUE)
      sp_max <- max(sp_vals, na.rm = TRUE)
    }
    
    sp_dt <- data.table(
        species = sp_name, 
        min = sp_min,
        mean = sp_mean,
        median = sp_median,
        q1 = sp_q1,
        q2 = sp_q2,
        q3 = sp_q3,
        max = sp_max, 
        filename = spec.filename
      )
    
    print(head(sp_vals, 1000))
    
    vebprint(sp_dt, veb = verbose, "species values:")
    
    return(sp_dt)
  }
  
  return(list(
    execute = execute
  ))
}

get_region_richness <- function(spec.filename, region, extra, verbose) {
  
  execute <- function(spec.filename, region, extra, verbose) {
    catn("Initiating data table.")
    sp_name <- gsub("-", " ", basename(dirname(spec.filename)))
    vebcat(sp_name, veb = verbose)
    
    sp_rast <- load_sp_rast(spec.filename)
    
    sp_region_dt <- extract_raster_to_dt(sp_rast, value = "cellAbundance", cells = TRUE)
    vebprint(sp_region_dt, veb = verbose)
    
    sp_region_dt <- sp_region_dt[!is.na(cellAbundance), ]
    
    sp_region_dt <- sp_region_dt[cellAbundance > 0]
    
    sp_region_dt[, cellAbundance := ifelse(any(cellAbundance > 0), 1, cellAbundance), by = "cell"]
    
    sp_region_dt[, species := sp_name]
    
    print(sp_region_dt)
    
    return(sp_region_dt)
  }
  
  process <- function(parallel.res, extra, verbose) { # parallel.res is a list of data tables
    res_list <- c(list(parallel.res$init.res), parallel.res$exec.res)
    
    # Read the cell csv made earlier
    cell_dt <- fread(extra)
    cell_dt[, richness := NULL]
    
    region_names <- list(
      subRegionCode = "FLOREG",
      subRegionName = "floregName",
      subRegionLong = "floregLong",
      subRegionLat = "floregLat"
    )
    
    for (i in seq_along(region_names)) {
      setnames(cell_dt, region_names[[i]], names(region_names)[i])
    }
    
    vebprint(cell_dt, veb = verbose, "Cell data table:")
    
    summed_list <- list()
    
    for (i in 1:length(res_list)) {
      cat("\rProcessing list", i, "/", length(res_list))
      
      sp_dt <- res_list[[i]]
      
      merged_dt <- merge(sp_dt, cell_dt, by = "cell")
      
      merged_dt[, regionAbundance := sum(cellAbundance), by = "subRegionCode"]
      
      merged_dt <- unique(merged_dt, by = "subRegionCode")
      
      summed_list[[i]] <- merged_dt
      
    }; catn()
    
    
    combined_dt <- rbindlist(summed_list)
    
    # Calculate region richness
    combined_dt[, regionRichness := uniqueN(species[cellAbundance > 0], na.rm = TRUE), by = "subRegionCode"]
    
    # Calculate country richness
    combined_dt[, countryRichness := uniqueN(species[cellAbundance > 0], na.rm = TRUE), by = "country"]
    
    return(combined_dt)
  }
  
  return(list(
    execute = execute,
    process = process
  ))
}

append_taxon <- function(dt, dt.append, dt.species, verbose = FALSE) {
  vebcat("Appending taxon names to data table.", color = "funInit")
  
  setnames(dt, "species", "cleanName")
  
  merged_dt <- merge(dt, dt.append, by = "cleanName")
  
  vebprint(dt, veb = verbose)
  vebcat("Number of species:", length(unique(dt[[dt.species]])), veb = verbose)
  
  species <- unique(dt[[dt.species]])
  
  species <- gsub("-", " ", species)
  
  catn("Checking GBIF backbone.\n")
  checklist <- name_backbone_checklist(name_data = species)
  
  checklist <- data.table(checklist[, c("kingdom", "phylum", "order", "family", "genus", "canonicalName", "verbatim_name")])
  
  names(checklist)[7] <- "species"
  
  na_sp <- checklist[is.na(checklist$genus),]
  
  if (nrow(na_sp) > 0) {
    vebcat("Species missing classification.", color = "nonFatalError")
    print(na_sp)
  }
  
  checklist[, species := gsub(" ", "-", species)]
  
  vebprint(checklist, veb = verbose)
  vebprint(dt, veb = verbose)
  
  catn("Merging data tables.")
  merged_dt <- merge(dt, checklist, by = "species")
  
  merged_dt <- data.table(merged_dt)
  
  vebcat("Taxon names successfully appended to data table.", color = "funSuccess")
  
  return(merged_dt)
}

calculate_taxon_richness <- function(dt, taxon, verbose = FALSE) {
  catn("Calculating richness for", highcat(taxon))
  
  taxon_col <- dt[[taxon]]
  
  dt[, taxonRichness := uniqueN(taxon_col[cellAbundance > 0], na.rm = TRUE), by = "subRegionCode"]
  
  # Calculate total taxon richness
  dt[, totalRichness := sum(taxonRichness, na.rm = TRUE), by = "subRegionCode"]
  
  # Calculate relative richness
  dt[, relativeRichness := taxonRichness / totalRichness]
  
  return(dt)
}

append_src_country <- function(src.dt, target.dt, level, taxon, verbose) {
  vebcat("Appending source regions to target data table", color = "funInit")
  catn("Getting working groups.")
  tdwg2 <- vect("./resources/region/tdwg/level2/level2.shp")
  tdwg2_dt <- convert_spatial_dt(tdwg2)
  tdwg2_dt <- tdwg2_dt[, .(LEVEL2_NAM, LEVEL2_COD)]
  colnames(tdwg2_dt) <- c("lvl2Name", "lvl2Code")
  
  tdwg3 <- vect("./resources/region/tdwg/level3/level3.shp")
  tdwg3_dt <- convert_spatial_dt(tdwg3)
  tdwg3_dt <- tdwg3_dt[, .(LEVEL3_NAM, LEVEL3_COD, LEVEL2_COD)]
  colnames(tdwg3_dt) <- c("lvl3Name", "lvl3Code", "lvl2Code")
  
  merged_tdwg <- merge(tdwg3_dt, tdwg2_dt, by = "lvl2Code")

  tdwg4 <- vect("./resources/region/tdwg/level4/level4.shp")
  tdwg4_dt <- convert_spatial_dt(tdwg4)
  tdwg4_dt <- tdwg4_dt[, .(Level_4_Na, Level3_cod)]
  colnames(tdwg4_dt) <- c("srcCountry", "lvl3Code")
  
  merged_tdwg <- merge(merged_tdwg, tdwg4_dt, by = "lvl3Code")

  catn("Setting up source country.")
  
  src_dt <- copy(src.dt)
  target_dt <- copy(target.dt)
  
  # Get source country
  src_dt <- src_dt[, .(species, srcCountry = country)]
  
  # combine with tdwg
  src_dt <- merge(src_dt, merged_tdwg, by = "srcCountry")
  
  if (taxon == "species") {
    src_dt <- unique(src_dt, by = "species")
  }

  print(src_dt)
  
  catn("Setting up target region.")
  
  target_dt[, taxonSum := uniqueN(.SD[[taxon]][value == 1], na.rm = TRUE), by = FLOREG]
  
  target_dt <- unique(target_dt, by = c(taxon, "FLOREG"))
  
  target_dt <- target_dt[taxonSum > 0]
  
  catn("Merging data tables.")
  sankey_dt <- merge(target_dt, src_dt, by = taxon, allow.cartesian = TRUE, all = TRUE)
  
  sankey_dt <- unique(sankey_dt, by = c("FLOREG", "srcCountry"))
  
  vebcat("Appended source regions successfully", color = "funSuccess")
  
  return(sankey_dt)
}

get_connections <- function(dt, subset, verbose = FALSE) {
  dt <- dt[, .(cleanName, order, family, genus, subRegionName, subRegionLong, subRegionLat, regionRichness, originCountry, originMeanLong, originMeanLat, taxonRichness, relativeRichness)]
  
  dt <- dt[!is.na(subRegionName) & subRegionName != "" & !is.na(originCountry) & originCountry != ""]
  
  dt <- unique(dt, by = subset)
  
  return(dt)
}

get_con_points <- function(dt, projection, longitude, latitude, verbose = FALSE) {
  
  vebprint(dt, verbose, "Input dt:")
  
  points <- vect(dt, geom=c(longitude, latitude), crs = longlat_crs)
  
  vebprint(points, verbose, "Points data:")
  
  points <- reproject_region(points, projection)
  
  vebprint(points, verbose, "Reprojected points:")

  return(points)
}
