clean_output_species <- function(out.dir, dt.result, dt.expected, clean.dirs = NULL, clean.files = NULL, verbose = FALSE) {
  cleaned_file <- paste0(out.dir, "/cleaned-results.csv")
  unexpected_file <- paste0(out.dir, "/unexpected-results.csv")
  
  if (file.exists(cleaned_file) && file.exists(unexpected_file)) {
    catn("Reading cleaned file.")
    cleaned_results <- fread(cleaned_file)
  } else {
    vebcat("Cleaning results.", color = "funInit")
    
    result_dt <- fread(dt.result)
    expected_dt <- fread(dt.expected, sep = "\t")
    
    result_names <- gsub("\\s+", "-", result_dt$cleanName)
    expected_names <- gsub("\\s+", "-", expected_dt$scientificName)
    
    # Clean the result file
    unmatched <- result_names[!result_names %in% expected_names]
    unmatched <- sort(unmatched)
    unmatched_check <- gsub("-", " ", unmatched) # cahnge - to _
    # Get rows that are not expected to exist in the data
    unexpected_results <- result_dt[result_dt$cleanName %in% unmatched_check, ]
    catn("Writing unexpected results to:", colcat(unexpected_file, color = "output"))
    fwrite(unexpected_results, unexpected_file, bom = TRUE)
    
    # Get rows that are expected to exist in the data
    cleaned_results <- result_dt[!result_dt$cleanName %in% unmatched_check, ]
    catn("Writing cleaned results to:", colcat(cleaned_file, color = "output"))
    fwrite(cleaned_results, cleaned_file, bom = TRUE)
    
    # Clean result directories
    if (is.null(clean.dirs) && is.null(clean.files)) {
      vebcat("Skipping cleaning of directories and files.", veb = verbose)
    } else {
      if (!is.null(clean.dirs)) {
        
        dirs <- list.dirs(clean.dirs)
        
        matched_dirs <- sapply(dirs, function(d) {
          base <- gsub("\\s+", "-", basename(d)) # change - to _
          
          if (base %in% unmatched) {
            return(d)
          } else {
            return(NA)
          }
        })
        
        matched_dirs <- matched_dirs[!is.na(matched_dirs)]
        matched_dirs <- unname(matched_dirs)
        
        vebprint(matched_dirs, verbose, "Matched Directories:")
        
        catn("Number of matched direcotries:", highcat(length(matched_dirs)))
        unlink(matched_dirs, recursive = TRUE)
        
        post_dirs <- list.dirs(clean.dirs)
        
        l_diff <- length(dirs) - length(post_dirs)
        
        if (l_diff == 0) {
          vebcat("Error, could not remove directories.", color = "nonFatalError")
        } else {
          vebcat("Successfully removed", l_diff, "directories.", color = "proSuccess")
        }
        
      } else {
        vebcat("Skipping cleaning of directories.", veb = verbose)
      }
      
      if (!is.null(clean.files)) {
        files <- list.files(clean.files)
        
        matched_files <- sapply(files, function(f) {
          base <- gsub("\\s+", "_", basename(f))
          
          if (base %in% unmatched) {
            return(f)
          } else {
            return(NA)
          }
        })
        
        matched_files <- matched_files[!is.na(matched_files)]
        matched_files <- unname(matched_files)
        
        vebprint(matched_files, verbose, "Matched Files:")
        
        catn("Number of matched Files:", highcat(length(matched_files)))
        file.remove(matched_files)
        
        post_files <- list.files(clean.files)
        
        l_diff <- length(files) - length(post_files)
        
        if (l_diff == 0) {
          vebcat("Error, could not remove files.", color = "nonFatalError")
        } else {
          vebcat("Successfully removed", l_diff, "files.", color = "proSuccess")
        }
        
      } else {
        vebcat("Skipping cleaning of files.", veb = verbose)
      }
    }

    vebcat("Results successfully cleaned.", color = "funSuccess")
    return(cleaned_results)
  }
}


filter_stats_data <- function(vis.dt, out.dir, hv.dir, hv.method, verbose = FALSE) {
  vebcat("Initializing stats data.", color = "funInit")
  
  inc_sp_out <- paste0(out.dir, "/included-species.csv")
  exc_sp_out <- paste0(out.dir, "/excluded-species.csv")
  
  if (!file.exists(inc_sp_out) || !file.exists(exc_sp_out) ) {
    
    included_sp <- vis.dt[excluded == FALSE, ]
    
    fwrite(included_sp, inc_sp_out, bom = TRUE)
    
    vebcat("Number of included species:", highcat(length(unique(included_sp$cleanName))), veb = verbose)
    
    excluded_sp <- vis.dt[excluded == TRUE, ]
    
    fwrite(excluded_sp, exc_sp_out, bom = TRUE)
    
    vebcat("Number of excluded species:", highcat(length(unique(excluded_sp$cleanName))), veb = verbose)
    
    mdwrite(
      post_seq_nums, 
      heading = "1;Visualize Sequence",
    )
    
    mdwrite(
      post_seq_nums, 
      heading = paste0("Number of included species: **", length(unique(included_sp$cleanName)), "**  \n",
                       "Number of Excluded species: **", length(unique(excluded_sp$cleanName)), "**"),
    )
    
  } else {
    excluded_sp <- fread(exc_sp_out)
    included_sp <- fread(inc_sp_out)
  }
  
  rm(excluded_sp)
  invisible(gc())
  
  vebcat("Stats successfully initialised.", color = "funSuccess")
  
  return(included_sp)
}

get_region_cells <- function(shape, template.filename, out.dir, verbose = FALSE) {
  
  vebcat("Converting region to data table with raster extents.", color = "funInit")
  
  out_file <- paste0(out.dir, "/region-cell.csv")
  
  if (file.exists(out_file)) {
    cell_regions_dt <- fread(out_file)
  } else {
    region <- load_region(shape)
    
    sp_rast <- load_sp_rast(template.filename) # Use as template
    
    if (!identical(ext(sp_rast), ext(region))) {
      catn("Cropping region extent to match species")
      region <- crop(region, ext(sp_rast))
    }
    
    catn("Extracting raster from region.")
    sp_cell_dt <- extract_raster_to_dt(sp_rast, region, value = "toRemove", cells = TRUE)
    sp_cell_dt[, toRemove := NULL]
    sp_cell_dt <- unique(sp_cell_dt, by = "cell")
    
    vebprint(sp_cell_dt, verbose)
    catn("Difference between nrows and unique cells:", highcat(nrow(sp_cell_dt) - length(unique(sp_cell_dt$cell))))
    
    catn("Converting region to data table.")
    region_dt <- as.data.frame(region)
    region_dt <- as.data.table(region_dt)
    region_dt[, ID := .I]
    region_dt <- handle_region_dt(region_dt)
    
    vebprint(region_dt, verbose)
    
    catn("Merging raster_dt and region_dt.")
    sp_regions_dt <- merge(sp_cell_dt, region_dt, by = "ID")
    setnames(sp_regions_dt, "ID", "regionId")
    
    sp_regions_dt <- unique(sp_regions_dt, by = "cell")
    
    catn("Difference between nrows and unique cells:", highcat(nrow(sp_regions_dt) - length(unique(sp_regions_dt$cell))))
    
    vebprint(sp_regions_dt, verbose)
    
    # Merge cell regions with the raster cells
    catn("Extracting raster cells.")
    sp_dt <- extract_raster_to_dt(sp_rast, value = "toRemove", cells = TRUE)
    sp_dt[, toRemove := NULL]
    
    vebprint(sp_dt, verbose)
    
    catn("Merging raster cells with region by cells.")
    cell_regions_dt <- merge(sp_dt, sp_regions_dt, by = "cell", all = FALSE)
    
    cell_regions_dt[country == "", country := NA]
    cell_regions_dt[subRegionName == "", subRegionName := NA]
    
    catn("Difference between nrows and unique cells:", highcat(nrow(cell_regions_dt) - length(unique(cell_regions_dt$cell))))
    
    vebprint(unique(cell_regions_dt$subRegionName), verbose, text = "SubRegions:")
    
    fwrite(cell_regions_dt, out_file, bom = TRUE)
  }
  
  vebcat("Final Number of subRegions:", highcat(length(unique(cell_regions_dt$subRegionName))))
  catn("Final Number of rows:", highcat(nrow(cell_regions_dt)))
  
  vebcat("Successfully converted region to data table with raster extents.", color = "funSuccess")
  
  return(cell_regions_dt)
}

get_inclusion_cell <- function(spec.filename, region = NULL, extra = NULL, verbose = FALSE) {
  
  init <- function(spec.filename, region, verbose) {
    sp_name <- basename(dirname(spec.filename[1]))
    vebcat("Species:", sp_name)
    
    sp_rast <- load_sp_rast(spec.filename[1])
    vebprint(sp_rast, text = "Raster:")
    
    summed_dt <- extract_raster_to_dt(sp_rast, value = "cellRichness", cells = TRUE)
    summed_dt <- unique(summed_dt, by = "cell")
    
    return(summed_dt)
  }
  
  execute <- function(spec.filename, region, extra = NULL, verbose) {
    
    # Do it in batches, that is faster and saves space
    sp_name <- basename(dirname(spec.filename[1]))
    vebcat("Species:", sp_name)
    
    sp_rast <- load_sp_rast(spec.filename[1])
    vebprint(sp_rast, text = "Raster:")
    
    summed_dt <- extract_raster_to_dt(sp_rast, value = "richnessSum", cells = TRUE)
    summed_dt <- unique(summed_dt, by = "cell")
    
    for (i in 2:length(spec.filename)) {
      cat("\rConverting raster to data table", i, "/", length(spec.filename), "| ")
      spec_file <- spec.filename[[i]]
      spec_name <- basename(dirname(spec.filename))
      
      vebcat("Species:", spec_name, veb = verbose)
      
      sp_rast <- load_sp_rast(spec_file)
      vebprint(sp_rast, verbose, text = "Raster:")
      
      res_dt <- extract_raster_to_dt(sp_rast, value = "cellAbundance", cells = TRUE)
      
      res_dt <- unique(res_dt, by = "cell")
      
      vebprint(res_dt, verbose, text = "Finished Raster:")
      
      summed_dt <- merge(summed_dt, res_dt, by = "cell")
      
      summed_dt[, richnessSum := ifelse(!is.na(richnessSum) & !is.na(cellAbundance), richnessSum + cellAbundance, ifelse(!is.na(richnessSum), richnessSum, cellAbundance)), by = cell]
      
      summed_dt[, cellAbundance := NULL]
      
      catn("nrow", nrow(summed_dt), "| max Richness", max(summed_dt$richnessSum, na.rm = TRUE))
    };catn()
    
    return(summed_dt)
  }
  
  # Handle the parallel results will return a list of results
  process <- function(parallel.res, extra, verbose) { # parallel.res is a list of data tables
    dt_summed <- parallel.res[[1]] 
    dt_list <- parallel.res[[2]]
    
    vebcat("Number of rows:", highcat(nrow(dt_summed)))
    vebcat("Number of unique cells:", highcat(length(unique(dt_summed$cell))))
    
    vebprint(dt_init, veb = verbose, "Init data table:")
    vebprint(dt_list, veb = verbose, "dt list:")
    
    for (i in 1:length(dt_list)) {
      cat("\rProcessing list", i, "/", length(dt_list))
      
      dt_sp <- dt_list[[i]] # cellRichness | richnessSum
      
      merged_dt <- merge(dt_summed, dt_sp, by = "cell")
      
      # Sum the richness
      merged_dt[, cellRichness := ifelse(!is.na(cellRichness) & !is.na(richnessSum), cellRichness + richnessSum, ifelse(!is.na(cellRichness), cellRichness, richnessSum)), by = cell]
      
      merged_dt[, richnessSum := NULL]
      
      dt_summed <- merged_dt
    }; catn()
    
    
    vebcat("Number of rows:", highcat(nrow(dt_summed)))
    vebcat("Number of unique cells:", highcat(length(unique(dt_summed$cell))))
    
    return(dt_summed)
  }
  
  return(list(
    init = init,
    execute = execute,
    process = process
  ))
}

convert_template_raster <- function(input.values, template.filename, projection, projection.method = "near", out.dir, verbose = F) {
  
  out_file <- paste0(out.dir, "/hotspots.tif")
  
  if (file.exists(out_file)) {
    catn("Reading existing raster.")
    template <- terra::rast(out_file)
    return(template)
  } 
  
  vebcat("Converting raster", color = "funInit")
  
  template <- terra::rast(template.filename)
  # Check length
  catn("Length of input / raster")
  catn(length(input.values), "/", ncell(template))
  
  terra::values(template) <- NA
  
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

get_paoo <- function(spec.filename, region, extra = NULL, verbose = FALSE) {
  
  execute <- function(spec.filename, region, extra, verbose) {
    sp_name <- basename(dirname(spec.filename))
    sp_name <- gsub("-", " ", sp_name)
    
    sp_rast <- terra::rast(spec.filename)
    
    catn("Converting region to data table.")
    region_dt <- as.data.frame(region)
    region_dt <- as.data.table(region_dt)
    region_dt[, ID := .I]
    region_dt <- handle_region_dt(region_dt)
    
    catn("Extracting raster cells.")
    sp_dt <- extract_raster_to_dt(sp_rast, region, value = "cellOccupancy", cells = TRUE)
    tot_cells <- sp_dt[!is.na(cell)]
    tot_cells <- uniqueN(tot_cells$cell)
    
    catn("Merging raster_dt and region_dt.")
    sp_regions_dt <- merge(sp_dt, region_dt, by = "ID")
    setnames(sp_regions_dt, "ID", "regionId")
    
    sp_regions_dt <- sp_regions_dt[!is.na(cellOccupancy)]
    
    cell_region_count <- sp_regions_dt[, nCellsRegion := uniqueN(cell, na.rm = TRUE), by = "subRegionName"]
    cell_region_count <- cell_region_count[, .(nCellsRegion, subRegionName)]
    cell_region_count <- unique(cell_region_count, by = "subRegionName")
    
    # subset to only get values above 0
    sp_regions_dt <- sp_regions_dt[cellOccupancy > 0]
    
    # Set all values to 1
    sp_regions_dt <- sp_regions_dt[, cellOccupancy := 1]
    
    sp_regions_dt <- unique(sp_regions_dt, by = "cell")
    
    sp_regions_dt <- sp_regions_dt[!is.na(cellOccupancy)]
    
    sp_regions_dt[country == "", country := NA]
    sp_regions_dt[subRegionName == "", subRegionName := NA]
    
    sp_regions_dt <- sp_regions_dt[, .(cell, cellOccupancy, subRegionName, country, subRegionLong, subRegionLat, westEast)]
    
    sp_regions_dt <- sp_regions_dt[, species := sp_name]
    
    sp_tpaoo <- sum(sp_regions_dt$cellOccupancy)
    
    sp_proptpaoo <- sp_tpaoo / tot_cells
    
    vebcat("Species Total Potential Area of Occupancy:", sp_tpaoo)
    
    sp_regions_dt <- sp_regions_dt[!is.na(subRegionName)]
    
    sp_tpaoo_region <- sp_regions_dt[, PAoO := sum(cellOccupancy, na.rm = TRUE), by = subRegionName]
    
    sp_tpaoo_region <- sp_tpaoo_region[, PAoO := ifelse(is.na(PAoO), 0, PAoO)]
    
    #sp_tpaoo_region[, propPAoO := PAoO / uniqueN(cell, na.rm = TRUE), by = "subRegionName"]
    
    sp_tpaoo_region <- merge(sp_tpaoo_region, cell_region_count, by = "subRegionName")
    
    sp_tpaoo_region <- sp_tpaoo_region[, propPAoO := PAoO / nCellsRegion, by = subRegionName]
    
    vebprint(sp_tpaoo_region, text = "Species Potential Area of Occupancy in each subregion:")
    out_dt <- data.table(species = sp_name, TPAoO = sp_tpaoo, propTPAoO = sp_proptpaoo, filename = spec.filename)
    
    out_dt <- merge(out_dt, sp_tpaoo_region, by = "species", all = TRUE)
    
    out_dt <- out_dt[, .(species, TPAoO, propTPAoO, PAoO, propPAoO, subRegionName, country, subRegionLong, subRegionLat, westEast, filename)]
    
    out_dt <- unique(out_dt, by = "subRegionName")
    
    vebprint(out_dt, text = "Finished table:")
    
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
  
  vebcat("Successfully stacked rasters", color = "funSuccess")
  
  return(stack)
}

get_prob_stats <- function(spec.filename, region, extra = NULL, verbose = FALSE) {
  
  execute <- function(spec.filename, region, extra, verbose) {
    sp_name <- basename(dirname(spec.filename))
    sp_name <- gsub("-", " ", sp_name)
    
    sp_rast <- terra::rast(spec.filename)
    
    catn("Converting region to data table.")
    region_dt <- as.data.frame(region)
    region_dt <- as.data.table(region_dt)
    region_dt[, ID := .I]
    region_dt <- handle_region_dt(region_dt)
    
    catn("Extracting raster cells.")
    sp_dt <- extract_raster_to_dt(sp_rast, region, value = "cellSuitability", cells = TRUE)
    
    sp_dt <- sp_dt[!is.na(cellSuitability)]
    
    catn("Merging raster_dt and region_dt.")
    sp_regions_dt <- merge(sp_dt, region_dt, by = "ID")
    setnames(sp_regions_dt, "ID", "regionId")
    sp_regions_dt <- unique(sp_regions_dt, by = "cell")
    sp_regions_dt <- sp_regions_dt[, .(cellSuitability, subRegionName, country)]
    sp_regions_dt <- sp_regions_dt[, species := sp_name]
    
    sp_Tmin <- min(sp_regions_dt$cellSuitability, na.rm = TRUE)
    sp_Tmean <- mean(sp_regions_dt$cellSuitability, na.rm = TRUE)
    sp_Tmedian <- median(sp_regions_dt$cellSuitability, na.rm = TRUE)
    # sp_Tq1 <- quantile(sp_regions_dt$cellSuitability, probs = 0.25, na.rm = TRUE)
    # sp_Tq2 <- quantile(sp_regions_dt$cellSuitability, probs = 0.5, na.rm = TRUE)
    # sp_Tq3 <- quantile(sp_regions_dt$cellSuitability, probs = 0.75, na.rm = TRUE)
    sp_Tmax <- max(sp_regions_dt$cellSuitability, na.rm = TRUE)
    
    vebcat("Species Total min Suitability:", sp_Tmin)
    vebcat("Species Total mean Suitability:", sp_Tmean)
    vebcat("Species Total median Suitability:", sp_Tmedian)
    # vebcat("Species Total 1. quantile Suitability:", sp_Tq1)
    # vebcat("Species Total 2. quantile Suitability:", sp_Tq2)
    # vebcat("Species Total 3. quantile Suitability:", sp_Tq3)
    vebcat("Species Total max Suitability:", sp_Tmax)
    
    total_dt <- data.table(
        species = sp_name, 
        totalMin = sp_Tmin,
        totalMean = sp_Tmean,
        totalMedian = sp_Tmedian,
        # q1 = sp_Tq1,
        # q2 = sp_Tq2,
        # q3 = sp_Tq3,
        totalMax = sp_Tmax, 
        filename = spec.filename
      )
    
    data_region <- sp_regions_dt[, regionMin := min(cellSuitability, na.rm = TRUE), by = subRegionName]
    data_region <- sp_regions_dt[, regionMean := mean(cellSuitability, na.rm = TRUE), by = subRegionName]
    data_region <- sp_regions_dt[, regionMedian := median(cellSuitability, na.rm = TRUE), by = subRegionName]
    # data_region <- sp_regions_dt[, regionq1 := quantile(cellSuitability, probs = 0.25, na.rm = TRUE), by = subRegionName]
    # data_region <- sp_regions_dt[, regionq2 := quantile(cellSuitability, probs = 0.5, na.rm = TRUE), by = subRegionName]
    # data_region <- sp_regions_dt[, regionq3 := quantile(cellSuitability, probs = 0.75, na.rm = TRUE), by = subRegionName]
    data_region <- sp_regions_dt[, regionMax := max(cellSuitability, na.rm = TRUE), by = subRegionName]
    
    data_region <- data_region[regionMax > 0, ]
    
    vebprint(data_region, text = "Species data for each subRegion:")
    
    out_dt <- merge(total_dt, data_region, by = "species", all = TRUE)
    
    out_dt <- out_dt[, .(species, totalMin, totalMean, totalMedian, totalMax, 
                         regionMin, regionMean, regionMedian, regionMax, 
                         subRegionName, country, filename)]
    
    out_dt <- unique(out_dt, by = "subRegionName")
    
    vebprint(out_dt, veb = verbose, "Finished data table:")
    
    rm(total_dt, data_region)
    invisible(gc())
    
    return(out_dt)
  }
  
  return(list(
    execute = execute
  ))
}

calculate_taxon_richness <- function(dt, taxon, verbose = FALSE) {
  catn("Calculating richness for", highcat(taxon))
  
  dt_copy <- copy(dt)
  
  dt_copy[, taxonRichness := uniqueN(cleanName, na.rm = TRUE), by = c(taxon, "subRegionName")]
  
  vebprint(dt_copy, verbose, text = "Taxon Richness before sum:")
  
  # Calculate total taxon richness
  dt_copy[, totalRichness := uniqueN(cleanName, na.rm = TRUE), by = subRegionName]
  
  # Calculate relative richness
  dt_copy[, relativeRichness := taxonRichness / totalRichness]
  
  vebprint(dt_copy, verbose, text = "Taxon Richness after sum:")
  
  return(dt_copy)
}

get_taxon_richness <- function(paoo.file, stats, taxon, verbose = FALSE) {
  vebcat("Getting region richness", color = "funInit")
  
  paoo_dt <- fread(paoo.file)
  sp_stats <- copy(stats)
  
  paoo_dt <- paoo_dt[, .(species, TPAoO, PAoO, subRegionName, country, subRegionLong, subRegionLat, westEast)]
  
  setnames(paoo_dt, "species", "cleanName")
  
  old_names <- c("country", "countryCode", "meanLong", "meanLat") # Later changed to median
  new_names <- c("originCountry", "originCountryCode", "originMeanLong", "originMeanLat") 
  setnames(sp_stats, old_names, new_names)
  sp_stats <- sp_stats[, .(cleanName, kingdom, phylum, class, order, family, genus, species, infraspecificEpithet, originCountryCode, originCountry, originMeanLong, originMeanLat)]
  
  merged_dt <- merge(paoo_dt, sp_stats, by = "cleanName", allow.cartesian = TRUE)
  
  merged_dt <- merged_dt[subRegionName == "", subRegionName := NA]
  
  vebcat("Number of species:", highcat(length(unique(merged_dt$species))))
  vebcat("Number of subRegions:", highcat(length(unique(merged_dt$subRegionName))))
  
  merged_dt <- merged_dt[!is.na(subRegionName)]
  
  vebcat("Number of species after removing NA subregions:", highcat(length(unique(merged_dt$species))))
  vebcat("Number of subRegions after removing NA subregions:", highcat(length(unique(merged_dt$subRegionName))))
  
  vebprint(merged_dt, verbose, text = "Merged Data Table:")
  
  richness_dt <- calculate_taxon_richness(
    merged_dt, 
    taxon = taxon,
  )
  
  richness_dt <- get_order_group(
    richness_dt
  )
  
  richness_dt <- unique(richness_dt, by = c(taxon, "subRegionName"))
  
  richness_dt <- richness_dt[, groupRelativeRichness := sum(relativeRichness), by = .(group, subRegionName)]
  
  
  vebprint(richness_dt[, .(sumofRR = sum(relativeRichness)), by = "subRegionName"], verbose, text = "final sum of RelativeRichness per region:")
  
  vebcat("Successfully acquired region richness", color = "funSuccess")
  
  return(richness_dt)
}

get_connections <- function(dt, taxon, verbose = FALSE) {
  
  dt_copy <- copy(dt)
  
  # subset_names <- c(taxon, "group", "subRegionName", "subRegionLong", "subRegionLat", "originCountry", "originMeanLong", "originMeanLat", "taxonRichness", "relativeRichness")
  # 
  # dt_copy <- dt_copy[, ..subset_names]
  
  dt_copy <- dt_copy[!is.na(subRegionName) & subRegionName != "" & !is.na(originCountry) & originCountry != ""]
  
  dt_copy <- dt_copy[, connections := .N, by = c(taxon, "subRegionName", "originCountry")]
  
  dt_copy <- unique(dt_copy, by = c(taxon, "subRegionName", "originCountry"))
  
  return(dt_copy)
}

get_con_points <- function(dt, projection, longitude, latitude, verbose = FALSE) {
  
  vebprint(dt, verbose, "Input dt:")
  
  points <- vect(dt, geom=c(longitude, latitude), crs = longlat_crs)
  
  vebprint(points, verbose, "Points data:")
  
  points <- reproject_region(points, projection)
  
  vebprint(points, verbose, "Reprojected points:")

  return(points)
}