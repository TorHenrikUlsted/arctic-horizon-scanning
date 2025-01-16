check_output_species <- function(out.dir, dt.result, dt.expected, hv.removed, clean.dirs = NULL, clean.files = NULL, verbose = FALSE) {
  cleaned_file <- paste0(out.dir, "/cleaned-results.csv")
  unexpected_file <- paste0(out.dir, "/gbif-changed-results.csv")

  if (file.exists(unexpected_file)) {
    catn("Reading cleaned file.")
    result_dt <- fread(dt.result)
  } else {
    vebcat("Cleaning results.", color = "funInit")

    result_dt <- fread(dt.result)
    hv_removed <- fread(hv.removed)
    expected_dt <- fread(dt.expected, sep = "\t")
    
    result_names <- gsub("\\s+", config$species$file_separator, result_dt$cleanName)
    hv_removed <- gsub("\\s+", config$species$file_separator, hv_removed$species)
    expected_names <- gsub("\\s+", config$species$file_separator, expected_dt$scientificName)
    
    vebprint(unique(result_names), verbose, text = "Result:")
    vebprint(unique(hv_removed), verbose, text = "Removed:")
    vebprint(unique(expected_names),verbose, text = "Expected:")
    
    # get unmatched and write out
    unmatched <- result_names[!result_names %in% expected_names & !hv_removed %in% expected_names]
    unmatched <- gsub(config$species$file_separator, " ", unmatched)
    unmatched <- data.table(scientificName = unique(unmatched))
    
    catn("Writing GBIF changed names to:", colcat(unexpected_file, color = "output"))
    fwrite(unmatched, unexpected_file, bom = TRUE)
    
    # Check if included but no overlap
    failed_sp_ret <- unique(result_dt[excluded == FALSE & overlapRegion == 0]$cleanName)
    no_data_return <- length(failed_sp_ret)
    
    if (no_data_return > 0) {
      vebcat(highcat(no_data_return), "Species returned with no overlap, setting to excluded", color = "warning")
      
      vebprint(failed_sp_ret, text = "Species with no data that generated projections:")
      
      stop("Need to solve above mentioned issues. Remove projections and set excluded to be: stats_dt[excluded := fifelse(overlapRegion == 0, TRUE, FALSE)]")
    }

    vebcat("Results successfully checked", color = "funSuccess")
    return(result_dt)
  }
}


filter_stats_data <- function(vis.dt, known.list, unknown.chunk.dir, out.dir, hv.dir, hv.method, verbose = FALSE) {
  vebcat("Initializing stats data.", color = "funInit")
  
  inc_sp_out <- paste0(out.dir, "/included-species.csv")
  exc_sp_out <- paste0(out.dir, "/excluded-species.csv")
  infra_sp_out <- paste0(out.dir, "/included-infraspecifics-region.csv")
  sp_w_infra_out <- paste0(out.dir, "/species-with-infraspecifics.csv")
  
  if (any(!sapply(c(inc_sp_out, exc_sp_out, infra_sp_out, sp_w_infra_out), file.exists))) {
    included_sp <- vis.dt[excluded == FALSE, ]
    
    catn("Writing out to file:", colcat(inc_sp_out, color = "output"))
    
    fwrite(included_sp, inc_sp_out, bom = TRUE)
    
    vebcat("Number of included species:", highcat(length(unique(included_sp$cleanName))))
    
    # Get underlying taxons and check if they are in the analysis
    
    infrasp <- copy(included_sp)
    
    setnames(infrasp, "cleanName", "species")
    
    infrasp <- unique(infrasp, by = "species")
    
    infrasp <- infrasp[, .(species)]
    
    known_sp <- fread(known.list, sep = "\t")
    
    catn("Acquiring underlying taxon ranks.")
    
    infrasp[, matched_infraSpecificEpithet := sapply(species, function(x) {
      # cat("\rProcessing species", x)
      matched_subspecies <- known_sp$scientificName[grepl(x, known_sp$scientificName)]
      if (length(matched_subspecies) > 0) matched_subspecies else NA
    })]
    
    infrasp <- infrasp[!is.na(matched_infraSpecificEpithet), ]
    
    infrasp <- infrasp[, list(matched_infraSpecificEpithet = unlist(matched_infraSpecificEpithet)), by = species]
    
    vebcat("Found", highcat(length(unique(infrasp$species))), "species with underlying taxons.")
    
    catn("Writing out to file:", colcat(infra_sp_out, color = "output"))
    fwrite(infrasp, infra_sp_out, bom = TRUE)
    
    # Check if they are in the analysis
    unknwon_files <- list.files(unknown.chunk.dir)
    
    unknwon_names <- gsub(config$species$file_separator, " ", gsub(".csv", "", unknwon_files))
    
    matches <- unique(infrasp$species) %in% unknwon_names
    
    matched <- unique(infrasp$species)[matches]
    
    infraspecifics <- data.table(species = character(), taxonRank = character(), infraspecificEpithet = character(), scientificName = character())
    
    if (length(matched > 0)) {
      for (i in 1:length(matched)) {
        cat("\rChecking species if included in analysis:", i, "/", length(matched))
        sp_file <- paste0(unknown.chunk.dir, "/", paste0(gsub(" ", config$species$file_separator, matched[i]), ".csv"))
        sp_name <- matched[i]
        sp <- fread(sp_file)
        
        infrasp <- unique(sp$infraspecificEpithet)
        scientificName <- unique(sp$scientificName)
        taxonRank <- unique(sp$taxonRank)
        
        infraspecifics <- rbind(infraspecifics, data.table(
          species = sp_name,
          taxonRank = taxonRank,
          infraspecificEpithet = infrasp,
          scientificName = scientificName
        ))
      }
      catn()
    }
    
    vebcat("Found", highcat(length(unique(!is.na(infraspecifics$infraspecificEpithet)))), "Species with additional underlying taxons.")
    
    catn("Writing out to file:", colcat(sp_w_infra_out, color = "output"))
    fwrite(infraspecifics, sp_w_infra_out, bom = TRUE)
    
    # Get excluded results
    
    excluded_sp <- vis.dt[excluded == TRUE, ]
    
    vebcat("Number of excluded species:", highcat(length(unique(excluded_sp$cleanName))))
    
    catn("Writing out to file:", colcat(exc_sp_out, color = "output"))
    fwrite(excluded_sp, exc_sp_out, bom = TRUE)
    
    mdwrite(
      config$files$post_seq_md,
      text = "1;Visualize Sequence",
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0(
        "Number of included species: **", length(unique(included_sp$cleanName)), "**  \n",
        "Number of Excluded species: **", length(unique(excluded_sp$cleanName)), "**"
      ),
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

split_spec_by_group <- function(spec, match.dt = NULL, match.colname = NULL, is.file = FALSE, verbose = FALSE) {
  
  if (!is.null(match.dt) && is.null(match.colname) || is.null(match.dt) && !is.null(match.colname)) {
    vebcat("match.dt and match.colname cannot have one as NULL when being used.", color = "fatalError")
    stop("Edit split_spec_by_group(match.dt, match.colname)")
  }
  
  if (is.character(spec)) {
    if (is.file) {
      spec <- data.table(filename = spec)
    } else {
      spec <- data.table(species = spec)
    }
  }
  
  match.dt <- unique(match.dt[, .(get(match.colname), order)], by = "V1")
  setnames(match.dt, "V1", "species")
  
  if (is.data.table(spec)) {
    init_length <- nrow(spec)
    if (is.file) spec[, species := clean_spec_filename(filename)]
    spec_taxons <- spec[match.dt, on = "species"]
  }
  
  vebprint(spec_taxons, verbose, "Species data table with order:")
  
  spec_group <- get_spec_group_dt(spec_taxons, "species")
  
  setcolorder(spec_group, setdiff(names(spec_group), "filename"))
  
  vebprint(spec_group, verbose, "Species data table with groups:")
  
  spec_group <- split(spec_group, by = "group")
  
  a_n <- if ("angiosperms" %in% names(spec_group)) nrow(spec_group$angiosperms) else 0
  g_n <- if ("gymnosperms" %in% names(spec_group)) nrow(spec_group$gymnosperms) else 0
  p_n <- if ("pteridophytes" %in% names(spec_group)) nrow(spec_group$pteridophytes) else 0
  
  res_length <- a_n + g_n + p_n
  
  if (init_length < res_length) {
    vebcat("Initial number of species", highcat(init_length), "is less than resulting length", highcat(res_length), color = "nonFatalError")
  }
  
  return(spec_group)
}

combine_groups <- function(x, out.order, out.n = NULL) {
  
  group <- rbindlist(x)
  
  if (grepl("-", out.order)) {
    out.order <- gsub("-", "", out.order)
    indecies <- order(-group[[out.order]])
    group <- group[indecies, ]
  } else {
    indecies <- order(group[[out.order]])
    group <- group[indecies, ]
  }
  
  if (!is.null(out.n)) group <- group[1:out.n]
  
  return(group)
}

get_inclusion_cell <- function(spec.filename, region = NULL, extra = NULL, verbose = FALSE) {
  init <- function(spec.filename, region, verbose) {
    sp_name <- clean_spec_filename(dirname(spec.filename[1]))
    vebcat("Species:", sp_name)
    
    sp_rast <- load_sp_rast(spec.filename[1])
    vebprint(sp_rast, text = "Raster:")
    
    summed_dt <- extract_raster_to_dt(sp_rast, value = "cellRichness", cells = TRUE)
    summed_dt <- unique(summed_dt, by = "cell")
    
    return(summed_dt)
  }
  
  execute <- function(spec.filename, region = NULL, extra = NULL, verbose = NULL) {
    # Do it in batches, that is faster and saves space
    sp_name <- clean_spec_filename(dirname(spec.filename[1]))
    vebcat("Species:", sp_name)
    
    sp_rast <- load_sp_rast(spec.filename[1])
    vebprint(sp_rast, text = "Raster:")
    
    summed_dt <- extract_raster_to_dt(sp_rast, value = "richnessSum", cells = TRUE)
    summed_dt <- unique(summed_dt, by = "cell")
    
    if (length(spec.filename) < 2) {
      return(summed_dt)
    }
    
    setkey(summed_dt, cell)
    
    for (i in 2:length(spec.filename)) {
      cat("\rConverting raster to data table", i, "|", length(spec.filename), "| ")
      spec_file <- spec.filename[[i]]
      spec_name <- clean_spec_filename(dirname(spec.filename))
      vebcat("Species:", spec_name, veb = verbose)
      
      sp_rast <- load_sp_rast(spec_file)
      vebprint(sp_rast, verbose, text = "Raster:")
      
      res_dt <- extract_raster_to_dt(sp_rast, value = "cellAbundance", cells = TRUE)
      res_dt <- unique(res_dt, by = "cell")
      
      vebprint(res_dt, verbose, "Finished Raster:")
      
      summed_dt[res_dt, on = "cell", richnessSum := 
                  fifelse(is.na(richnessSum) & is.na(i.cellAbundance), NA_real_,
                          fifelse(is.na(richnessSum), i.cellAbundance,
                                  fifelse(is.na(i.cellAbundance), richnessSum,
                                          richnessSum + i.cellAbundance)))
      ]
      
      catn("nrow", nrow(summed_dt), "| max Richness", max(summed_dt$richnessSum, na.rm = TRUE))
    }
    catn()
    
    return(summed_dt)
  }
  
  # Handle the parallel results will return a list of results
  process <- function(parallel.res, extra, verbose) { # parallel.res is a list of data tables
    dt_summed <- parallel.res[[1]]
    dt_list <- parallel.res[[2]]
    
    vebcat("Number of rows:", highcat(nrow(dt_summed)))
    vebcat("Number of unique cells:", highcat(length(unique(dt_summed$cell))))
    
    vebprint(dt_summed, veb = verbose, "Init data table:")
    vebprint(dt_list, veb = verbose, "dt list:")
    
    for (i in 1:length(dt_list)) {
      cat("\rProcessing list", i, "/", length(dt_list))
      
      dt_sp <- dt_list[[i]] 
      
      dt_summed[dt_sp, on = "cell", cellRichness := 
                  fifelse(is.na(cellRichness) & is.na(i.richnessSum), NA_real_, # NA + NA = NA
                          fifelse(is.na(cellRichness), i.richnessSum, # NA + value = value
                                  fifelse(is.na(i.richnessSum), cellRichness, # value + NA = value
                                          cellRichness + i.richnessSum))) # value1 + value2 = sum of values
      ]
      
    }
    catn()
    
    
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

convert_template_raster <- function(input.values, template.filename, projection, projection.method = "near", out.dir, verbose = FALSE) {
  if (grepl("/", out.dir)) {
    out_file <- paste0(out.dir, "-hotspots.tif")
  } else {
    out_file <- paste0(out.dir, "/hotspots.tif")
  }
  
  
  if (file.exists(out_file)) {
    catn("Reading", highcat(basename(out_file)))
    return(terra::rast(out_file))
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

get_world_map <- function(projection, map.type = "wgsrpd", scale = "level3", pole = NULL, verbose = FALSE) {
  vebcat("Setting up world Map.", color = "funInit")
  
  accepted_maps <- c("wgsrpd", "rnaturalearth")
  
  if (!(tolower(map.type) %in% accepted_maps)) {
    vebcat("Error: map.type not accepted.", color = "fatalError")
    catn("Accepted map.type:", accepted_maps)
    stop("Change map.type paramter.")
  }
  
  if (tolower(map.type) == "wgsrpd") {
    if (!grepl("level", scale)) {
      vebcat("Error: scale has to be in the form paste0('level', number), where number has to be either 1 (continents), 2 (regions), 3 (botanical countries), 4 (Basic Recording Units)", color = "fatalError")
      catn("For help see:", highcat("http://www.tdwg.org/standards/109"))
    }
    path <- paste0("./resources/region/wgsrpd/", scale)
    download_github_dir_if("tdwg", "wgsrpd", "master", scale, path)
    world_map <- load_region(paste0(path, "/", scale, ".shp"))
  } else if (tolower(map.type) == "rnaturalearth") {
    if (grepl("level", scale)) {
      vebcat("Error: scale has to be in the form 'large', 'medium', or 'small'.", color = "fatalError")
      catn("For help type: ?rnaturalearth")
    }
    world_map <- ne_countries(scale = scale, returnclass = "sf")
    world_map <- vect(world_map)
  } 
  
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
  
  print(projection)
  
  world_map <- suppressWarnings(check_crs(world_map, projection, verbose))
  
  vebcat("World Map setup completed successfully", color = "funSuccess")
  
  return(world_map)
}

get_paoo <- function(spec.filename, region, extra = NULL, verbose = FALSE) {
  execute <- function(spec.filename, region, extra, verbose) {
    sp_name <- basename(dirname(spec.filename))
    sp_name <- gsub(config$species$file_separator, " ", sp_name)
    
    sp_rast <- terra::rast(spec.filename)
    
    catn("Converting region to data table.")
    region_dt <- as.data.frame(region)
    region_dt <- as.data.table(region_dt)
    region_dt[, ID := .I]
    region_dt <- handle_region_dt(region_dt)
    
    catn("Extracting raster cells.")
    sp_dt <- extract_raster_to_dt(sp_rast, region, value = "cellOccupancy", cells = TRUE)
    tot_cells <- sp_dt[!is.na(cell), uniqueN(cell)]
    
    catn("Merging raster_dt and region_dt.")
    sp_regions_dt <- sp_dt[region_dt, on = "ID", nomatch = 0]
    #sp_regions_dt <- merge(sp_dt, region_dt, by = "ID")
    setnames(sp_regions_dt, "ID", "regionId")
    
    sp_regions_dt[, `:=`(
      cellOccupancy = fifelse(cellOccupancy > 0, 1, fifelse(cellOccupancy == 0, 0, NA_real_)),
      country = fifelse(country == "", NA_character_, country),
      subRegionName = fifelse(subRegionName == "", NA_character_, subRegionName),
      species = sp_name
    )]
    
    vebprint(any(!is.na(sp_regions_dt$cellOccupancy)), verbose, text = "Any not NA:")
    vebprint(any(sp_regions_dt$cellOccupancy != 0), verbose, text = "Any not 0:")
    
    sp_regions_dt <- unique(sp_regions_dt[!is.na(cellOccupancy) & !is.na(subRegionName)], by = "cell")
    
    sp_regions_dt[, nCellsRegion := uniqueN(cell, na.rm = TRUE), by = "subRegionName"]
    
    vebprint(sp_regions_dt, verbose, "cell_region_count:")
    
    sp_tpaoo <- sp_regions_dt[, sum(cellOccupancy)]
    sp_proptpaoo <- sp_tpaoo / tot_cells
    
    vebcat("Species Total Potential Area of Occupancy:", highcat(sp_tpaoo))
    
    sp_tpaoo_region <- sp_regions_dt[, .(
      PAoO = sum(cellOccupancy, na.rm = TRUE),
      nCellsRegion = .N,
      country = first(country),  # Assuming country is the same for each subRegionName
      species = first(species),   # Assuming species is the same for all rows
      subRegionLong = first(subRegionLong),
      subRegionLat = first(subRegionLat),
      westEast = first(westEast)
    ), by = "subRegionName"]
    
    sp_tpaoo_region[, propPAoO := PAoO / nCellsRegion, by = "subRegionName"]
    
    vebprint(sp_tpaoo_region, verbose, "Species Potential Area of Occupancy in each subregion:")
    
    out_dt <- sp_tpaoo_region[, .(
      species,
      TPAoO = sp_tpaoo,
      propTPAoO = sp_proptpaoo,
      PAoO,
      propPAoO,
      subRegionName,
      country,
      subRegionLong,
      subRegionLat,
      westEast,
      filename = spec.filename
    )]
    
    vebprint(out_dt, text = "Finished table:")
    
    return(out_dt)
  }
  
  return(list(
    execute = execute
  ))
}

stack_projections <- function(filenames, projection, projection.method, out.dir, binary = FALSE, verbose = FALSE) {
  vebcat("Stacking rasters", color = "funInit")
  
  method_dir <- file.path(out.dir, gsub(".tif", "", basename(filenames)[1]))
  
  nodes_dir <- file.path(method_dir, "nodes")
  rasters_dir <- file.path(method_dir, "rasters")
  
  create_dir_if(nodes_dir, rasters_dir)
  
  rast_out_names <- paste0(basename(dirname(filenames)), ".tif")
  
  raster_files <- list.files(rasters_dir, full.names = FALSE)
  
  skip <- FALSE
  
  if (all(rast_out_names %in% raster_files)) {
    catn("Found raster files, skipping to stacking.")
    skip <- TRUE
  }
  
  if (!skip) {
    
    tryCatch({
      max_cores <- calc_num_cores(
        ram.high = 5,
        verbose = FALSE
      )
      
      catn("Found", length(filenames), "files.")
      
      max_cores <- min(length(filenames), max_cores$high)
      
      catn("Creating cluster of", max_cores, "core(s).")
      
      cl <- makeCluster(max_cores)
      
      clusterEvalQ(cl, {
        source("./src/utils/utils.R")
        load_utils(parallel = TRUE)
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
      
      catn("Stacking rasters in parallel.")
      
      # Apply the function in parallel
      clusterApplyLB(cl, 1:length(filenames), function(i) {
        tryCatch({
          node_file <- paste0(nodes_dir, "/", paste0("node", i))
          create_file_if(node_file)
          
          try(node_con <- file(node_file, "a"))
          sink(node_con, type = "output")
          sink(node_con, type = "message")
          catn("Setting up", paste0("node-", i))
          
          projection <- get_crs_config(projection)
          
          file <- filenames_dt$filename[i]
          rast_out_names <- paste0(basename(dirname(file)), ".tif")
          name <- gsub(config$species$file_separator, " ", basename(dirname(file)))
          
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
          
          if (!identical(crs(raster, proj = TRUE), crs(projection, proj = TRUE))) {
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
        },
        error = function(e) {
          vebprint(e$message, text = "\nError message:")
          stop(e)
          
        }, finally = {
          catn("Parallel finished, cleaning up sink connections")
          sink(type = "message", NULL)
          sink(type = "output", NULL)
          if (exists("node_con") && isOpen(node_con)) {
            close(node_con)
          }
        })
      })
    }, error = function(e) {
      vebcat("An error occurred in the parallel process ~ stopping cluster and closing connections.", color = "fatalError")
      catn("Error messages written to node logs in:\n", colcat(nodes_dir, color = "output"))
      stop(e)
    }, finally =  {
      catn("Cleaning up cluster connections")
      if (exists("cl")) {
        tryCatch({
          stopCluster(cl)
        }, error = function(e) {
          vebcat("Error stopping cluster:", e$message, color = "fatalError")
        })
      }
      closeAllConnections()
    }) # clusterApplyLB
    
  } # skip
  
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
  }
  catn()
  
  vebcat("Successfully stacked rasters", color = "funSuccess")
  
  return(stack)
}

get_prob_stats <- function(spec.filename, region, extra = NULL, verbose = FALSE) {
  execute <- function(spec.filename, region, extra, verbose) {
    sp_name <- basename(dirname(spec.filename))
    sp_name <- gsub(config$species$file_separator, " ", sp_name)
    
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
    
    out_dt <- out_dt[, .(
      species, totalMin, totalMean, totalMedian, totalMax,
      regionMin, regionMean, regionMedian, regionMax,
      subRegionName, country, filename
    )]
    
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

calculate_taxon_richness <- function(dt, taxon, country = NULL, verbose = FALSE) {
  catn("Calculating richness for", highcat(taxon))
  
  dt_copy <- copy(dt)
  # Calculate richness per taxon and subregion (number of unique species in each combination)
  dt_copy[, taxonRichness := length(unique(cleanName[PAoO > 0])), 
          by = c(taxon, "subRegionName")]
  
  vebprint(dt_copy, verbose, text = "Taxon Richness before sum:")
  
  # Calculate total richness per subregion (across all taxa)
  dt_copy[, totalRichness := length(unique(cleanName[PAoO > 0])), 
          by = subRegionName]
  
  # Calculate relative richness
  dt_copy[, relativeRichness := taxonRichness / totalRichness, by = c(taxon, "subRegionName")]
  
  # Remove duplicates to get one row per taxon-subregion combination
  if (is.null(country)) {
    dt_copy <- unique(dt_copy, by = c(taxon, "subRegionName"))
  } else {
    dt_copy <- unique(dt_copy, by = c(taxon, "subRegionName", country))
  }
  
  vebprint(dt_copy, verbose, text = "Taxon Richness after sum:")
  
  return(dt_copy)
}

get_taxon_richness <- function(paoo.file, stats, taxon, country = FALSE, verbose = FALSE) {
  vebcat("Getting region richness", color = "funInit")
  
  if (taxon == "species") taxon = "cleanName"
  
  if (is.character(paoo.file)) {
    paoo.file <- fread(paoo.file)
  }
  sp_stats <- copy(stats)
  
  paoo_dt <- paoo.file[, .(species, TPAoO, PAoO, subRegionName, country, subRegionLong, subRegionLat, westEast)]
  
  setnames(paoo_dt, "species", "cleanName")
  
  old_names <- c("level3Name", "level3Code", "level3Long", "level3Lat")
  new_names <- c("originCountry", "originCountryCode", "originLong", "originLat")
  setnames(sp_stats, old_names, new_names)
  sp_stats <- sp_stats[, .(cleanName, kingdom, phylum, class, order, family, genus, species, infraspecificEpithet, originCountryCode, originCountry, originLong, originLat, level2Code, level1Code)]
  
  merged_dt <- paoo_dt[sp_stats, on = "cleanName", allow.cartesian = TRUE]
  
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
    country = if (country) "originCountry" else NULL,
    verbose = verbose
  )
  
  richness_dt <- get_order_group(
    richness_dt
  )
  
  if (!country) {
    richness_dt <- unique(richness_dt, by = c(taxon, "subRegionName"))
  } else {
    richness_dt <- unique(richness_dt, by = c(taxon, "subRegionName", "originCountry"))
  }
  
  richness_dt <- richness_dt[, groupRelativeRichness := sum(relativeRichness), by = .(group, subRegionName)]
  
  vebprint(richness_dt[, .(sumofRR = sum(relativeRichness)), by = "subRegionName"], verbose, text = "final sum of RelativeRichness per region:")
  
  vebcat("Successfully acquired region richness", color = "funSuccess")
  
  return(richness_dt)
}

get_connections <- function(dt, verbose = FALSE) {
  dt_copy <- copy(dt)
  
  dt_copy <- dt_copy[!is.na(subRegionName) & subRegionName != "" & !is.na(originCountry) & originCountry != ""]
  
  dt_copy <- dt_copy[, connections := data.table::uniqueN(cleanName), by = c("cleanName", "subRegionName", "originCountry")]
  
  dt_copy <- unique(dt_copy, by = c("cleanName", "subRegionName", "originCountry"))
  
  return(dt_copy)
}

get_con_points <- function(dt, projection, longitude, latitude, verbose = FALSE) {
  vebprint(dt, verbose, "Input dt:")
  
  points <- vect(dt, geom = c(longitude, latitude), crs = config$projection$crs$longlat)
  
  vebprint(points, verbose, "Points data:")
  
  points <- check_crs(points, projection, verbose = verbose)
  
  vebprint(points, verbose, "Reprojected points:")
  
  return(points)
}

convert_con_points <- function(dt, taxon, projection, centroid = FALSE, verbose = FALSE) {
  sub_dt <- copy(dt)
  
  if (taxon == "cleanName") {
    origin_subset <- sub_dt[, .(cleanName, originLong, originLat, connections, originCountry, level2Code)]
  } else {
    origin_subset <- sub_dt[, .(cleanName, get(taxon), originLong, originLat, connections, originCountry, level2Code)]
    setnames(origin_subset, "V2", taxon)
  }
  origin_subset <- unique(origin_subset, by = c("cleanName", "originCountry"))
  
  if (centroid) {
    origin_subset[, level3Name := originCountry]
    wgsrpd_centroids <- get_wgsrpd_polygon_centroid(origin_subset)
    origin_subset <- origin_subset[wgsrpd_centroids, on = "level3Name"]
    data.table::setnames(origin_subset, c("originLong", "originLat"), c("level3Long", "level3Lat"))
    data.table::setnames(origin_subset, c("centroidLong", "centroidLat"), c("originLong", "originLat"))
  }
  
  origin_points <- get_con_points(origin_subset, projection, "originLong", "originLat", verbose = verbose)
  
  if (taxon == "cleanName") {
    dest_subset <- sub_dt[, .(cleanName, subRegionLong, subRegionLat, subRegionName)]
  } else {
    dest_subset <- sub_dt[, .(cleanName, get(taxon), subRegionLong, subRegionLat, subRegionName)]
    setnames(dest_subset, "V2", taxon)
  }
  dest_subset <- unique(dest_subset, by = c("cleanName", "subRegionName"))
  dest_points <- get_con_points(dest_subset, projection, "subRegionLong", "subRegionLat", verbose = verbose)
  
  # Remove NA
  origin_points <- origin_points[!is.na(origin_points)]
  dest_points <- dest_points[!is.na(dest_points)]
  
  catn("Converting points back to data tables.")
  origin <- terra::as.data.frame(origin_points, geom = "xy")
  origin <- as.data.table(origin)
  setnames(origin, "x", "originX")
  setnames(origin, "y", "originY")
  origin <- unique(origin, by = c("cleanName", "originX", "originY"))
  
  dest <- terra::as.data.frame(dest_points, geom = "xy")
  dest <- as.data.table(dest)
  setnames(dest, "x", "destX")
  setnames(dest, "y", "destY")
  dest <- unique(dest, by = c("cleanName", "destX", "destY"))
  
  origin <- unique(origin, by = c("originX", "originY"))
  dest <- unique(dest, by = c("destX", "destY"))
  
  return(list(
    origin = origin,
    dest = dest
  ))
}

create_zoom_annotation <- function(basemap, region, points, verbose = FALSE) {
  # Get the bounding box of the subregion for the zoom
  bbox <- terra::ext(region)
  # Add some padding around the bbox
  padding <- c(
    (bbox[2] - bbox[1]) * 0.1,  # 10% padding
    (bbox[4] - bbox[3]) * 0.1
  )
  
  bbox <- c(
    bbox[1] - padding[1],
    bbox[2] + padding[1],
    bbox[3] - padding[2],
    bbox[4] + padding[2]
  )
  
  zoom_plot <- ggplot() +
    geom_spatvector(data = basemap) +
    geom_spatvector(
      data = region,
      fill = "#0013ff",
      color = "black",
      linewidth = 0.1
    ) +
    geom_point(
      data = points, 
      aes(x = originX, y = originY), 
      color = "darkgreen",
      size = 1, 
      stroke = 0, 
      shape = 16
    ) +
    coord_sf(
      xlim = c(bbox[1], bbox[2]),
      ylim = c(bbox[3], bbox[4])
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.background = element_rect(fill = "white", color = "black"),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  
  # Convert zoom plot to grob
  return(list(
    grob = ggplotGrob(zoom_plot),
    extents = bbox
  ))
}

position_zoom_annotation <- function(top, right, plot.ext, region.ext, width = NULL, height = NULL, verbose = FALSE) {
  # Calculate aspect ratio from the zoom region's extents
  zoom_x_range <- unname(region.ext[2]) - unname(region.ext[1])
  zoom_y_range <- unname(region.ext[4]) - unname(region.ext[3])
  aspect_ratio <- zoom_x_range / zoom_y_range
  
  # Get the plot ranges for positioning
  plot_x_range <- unname(plot.ext[2]) - unname(plot.ext[1])
  plot_y_range <- unname(plot.ext[4]) - unname(plot.ext[3])
  
  if (is.null(width)) width <- 0.25
  if (is.null(height)) height <- 0.25
  
  # Adjust dimensions to maintain aspect ratio while keeping the larger dimension fixed
  if (aspect_ratio > 1) {  # wider than tall
    base_width <- plot_x_range * width
    base_height <- plot_y_range * height
    actual_width <- base_width
    actual_height <- base_width / aspect_ratio
  } else {  # taller than wide
    base_width <- plot_x_range * (width + 0.1)
    base_height <- plot_y_range * (height + 0.1)
    
    actual_width <- base_height * aspect_ratio
    actual_height <- base_height
  }
  
  # Position consistently from top edge
  inset_ymax <- unname(plot.ext[4]) - (top * actual_height)
  inset_ymin <- inset_ymax - base_height  # Use base_height for consistent top spacing
  
  # Adjust ymin to maintain aspect ratio
  if (aspect_ratio > 1) {
    inset_ymin <- inset_ymax - actual_height
  }
  
  # Position from right edge
  inset_xmax <- unname(plot.ext[2]) - (right * plot_x_range)
  inset_xmin <- inset_xmax - actual_width
  
  list(
    xmin = inset_xmin,
    xmax = inset_xmax,
    ymin = inset_ymin,
    ymax = inset_ymax,
    aspect_ratio = aspect_ratio
  )
}

calc_list_rows <- function(dt.filenames, begin, end, write.md = FALSE) {
  mdwrite(
    config$files$post_seq_md,
    text = "1;List number of rows",
    veb = write.md
  )
  
  catn("Number of rows in the different files:")
  catn("--------------------------------------")
  
  for (i in 1:length(dt.filenames)) {
    file <- dt.filenames[i]
    
    name_dir <- progressive_dirname(file, begin, end)
    
    name <- paste0(basename(name_dir), " - ", basename(file))
    
    file_length <- system_calc_rows(file)
    
    catn(paste0(name, ":"), highcat(file_length))
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0("2;", name),
      veb = write.md
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0("Number of species: **", file_length, "**"),
      veb = write.md
    )
  }
  
  return(invisible())
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

split_spec_by_group <- function(spec, match.dt = NULL, match.colname = NULL, is.file = FALSE, verbose = FALSE) {
  
  if (!is.null(match.dt) && is.null(match.colname) || is.null(match.dt) && !is.null(match.colname)) {
    vebcat("match.dt and match.colname cannot have one as NULL when being used.", color = "fatalError")
    stop("Edit split_spec_by_group(match.dt, match.colname)")
  }
  
  if (is.character(spec)) {
    if (is.file) {
      spec <- data.table(filename = spec)
    } else {
      spec <- data.table(species = spec)
    }
  }
  
  match.dt <- unique(match.dt[, .(get(match.colname), order)], by = "V1")
  setnames(match.dt, "V1", "species")
  
  if (is.data.table(spec)) {
    init_length <- nrow(spec)
    if (is.file) spec[, species := clean_spec_filename(filename)]
    spec_taxons <- spec[match.dt, on = "species"]
  }
  
  vebprint(spec_taxons, verbose, "Species data table with order:")
  
  spec_group <- get_spec_group_dt(spec_taxons, "species")
  
  setcolorder(spec_group, setdiff(names(spec_group), "filename"))
  
  vebprint(spec_group, verbose, "Species data table with groups:")
  
  spec_group <- split(spec_group, by = "group")
  
  a_n <- if ("angiosperms" %in% names(spec_group)) nrow(spec_group$angiosperms) else 0
  g_n <- if ("gymnosperms" %in% names(spec_group)) nrow(spec_group$gymnosperms) else 0
  p_n <- if ("pteridophytes" %in% names(spec_group)) nrow(spec_group$pteridophytes) else 0
  
  res_length <- a_n + g_n + p_n
  
  if (init_length < res_length) {
    vebcat("Initial number of species", highcat(init_length), "is less than resulting length", highcat(res_length), color = "nonFatalError")
  }
  
  return(spec_group)
}

combine_groups <- function(x, out.order, out.n = NULL) {
  
  group <- rbindlist(x)
  
  if (grepl("-", out.order)) {
    out.order <- gsub("-", "", out.order)
    indecies <- order(-group[[out.order]])
    group <- group[indecies, ]
  } else {
    indecies <- order(group[[out.order]])
    group <- group[indecies, ]
  }
  
  if (!is.null(out.n)) group <- group[1:out.n]
  
  return(group)
}

latitude_analysis <- function(spec.filename, region = NULL, extra = NULL, verbose = FALSE) {
  
  execute <- function(spec.filename, region = NULL, extra = NULL, verbose = NULL) {
    sp_name <- clean_spec_filename(spec.filename)
    vebcat("Species:", sp_name)
    
    cols_to_sel <- c("cleanName", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", "countryCode", "stateProvince", "year")
    spec_dt <- fread(spec.filename, select = cols_to_sel)
    
    cleaned_data <- clean_spec_occ(
      dt = spec_dt,
      column = "cleanName", 
      projection = extra$projection, 
      resolution = config$projection$raster_scale_m, 
      seed = config$simulation$seed, 
      long = "decimalLongitude", 
      lat = "decimalLatitude",  
      verbose = verbose
    )
    
    if (is.null(cleaned_data)) {
      cleaned_data <- data.table(
        cleanName = sp_name,
        medianLat = NA_real_,
        nobs = 0,
        origNobs = nrow(spec_dt)
      )
      
      return(cleaned_data)
    }
    
    nobs <- nrow(cleaned_data)
    
    cleaned_data <- unique(cleaned_data, by = "cleanName")
    
    cleaned_data <- cleaned_data[, .(cleanName, medianLat)]
    
    cleaned_data[, `:=` (
      nobs = nobs, 
      origNobs = nrow(spec_dt)
    )]
    
    return(cleaned_data)
  }
  
  return(list(
    execute = execute
  ))
}

analyze_latitudinal_pattern <- function(stats, spec.file, out.dir, verbose = FALSE) {
  vebcat("Analyzing latitudinal patterns", color = "funInit")
  
  log_dir <- file.path(out.dir, "logs")
  out_file <- file.path(out.dir, "results/calculated-median-latitude.csv")
  lat_stats <- fread(stats)
  
  if (file.exists(out_file)) {
    catn("File found, reading file..")
    calced_dt <- fread(out_file)
  } else {
    catn("Getting species filenames")
    spec_name_stats <- unique(lat_stats$cleanName)
    
    spec_list <- readLines(spec.file)
    
    spec_names <- gsub(config$species$file_separator, " ", gsub(".csv", "", basename(spec_list)))
    
    spec_list <- spec_list[spec_names %in% spec_name_stats]
    
    projection <- get_crs_config("longlat")
    
    calced_dt <- parallel_spec_handler(
      spec.dirs = spec_list,
      dir = file.path(log_dir, "distribution"),
      hv.project.method = "longlat",
      out.order = "-medianLat",
      extra = list(projection = projection),
      fun = latitude_analysis,
      verbose = verbose
    )
  }
  
  lat_stats <- lat_stats[, .(cleanName, order, realizedNiche, overlapRegion, observations, dimensions, excluded)]
  
  lat_stats <- unique(lat_stats, by = "cleanName")
 
  catn("Species min absolute median latitude:", highcat(min(calced_dt$medianLat, na.rm = TRUE)))
  
  catn("Species max absolute median latitude:", highcat(max(calced_dt$medianLat, na.rm = TRUE)))
  
  merged_dt <- lat_stats[calced_dt, on = "cleanName", nomatch = 0L]
  
  merged_dt <- unique(merged_dt, by = "cleanName")
  
  merged_dt <- merged_dt[!is.na(medianLat)]
  
  return(merged_dt)
}

compare_gamlss_models <- function(data, predictor = "medianLat", response = "overlapRegion") {
  # Create formula strings
  full_formula <- as.formula(paste(response, "~", predictor))
  null_formula <- as.formula(paste(response, "~ 1"))
  
  # Start CLI output
  cli_h1("GAMLSS Model Comparison")
  cli_alert_info("Fitting models...")
  
  # List to store all models
  models <- list()
  
  # Progress bar
  cli_progress_bar("Fitting models", total = 5)
  
  # 1. Full model (latitude coefficient for mu, sigma, and nu)
  cli_progress_update()
  models$full <- gamlss(
    formula = full_formula,
    sigma.formula = as.formula(paste("~", predictor)),
    nu.formula = as.formula(paste("~", predictor)),
    family = BEZI,
    data = data,
    trace = FALSE
  )
  
  # 2. Latitude-zero only (latitude coefficient for nu, intercept for mu)
  cli_progress_update()
  models$zero_only <- gamlss(
    formula = null_formula,
    sigma.formula = as.formula(paste("~", predictor)),
    nu.formula = as.formula(paste("~", predictor)),
    family = BEZI,
    data = data,
    trace = FALSE
  )
  
  # 3. Latitude-magnitude only (latitude coefficient for mu, intercept for nu)
  cli_progress_update()
  models$magnitude_only <- gamlss(
    formula = full_formula,
    sigma.formula = as.formula(paste("~", predictor)),
    nu.formula = ~ 1,
    family = BEZI,
    data = data,
    trace = FALSE
  )
  
  # 4. Intercept only model (intercepts for both nu and mu)
  cli_progress_update()
  models$intercept_only <- gamlss(
    formula = null_formula,
    sigma.formula = as.formula(paste("~", predictor)),
    nu.formula = ~ 1,
    family = BEZI,
    data = data,
    trace = FALSE
  )
  
  # 5. Complete null model (intercepts for mu, sigma, and nu)
  cli_progress_update()
  models$complete_null <- gamlss(
    formula = null_formula,
    sigma.formula = ~ 1,
    nu.formula = ~ 1,
    family = BEZI,
    data = data,
    trace = FALSE
  )
  
  cli_progress_done()
  
  # Calculate AIC, BIC and deltas
  aic_values <- sapply(models, AIC)
  bic_values <- sapply(models, function(m) BIC(m))
  delta_aic <- aic_values - min(aic_values)
  delta_bic <- bic_values - min(bic_values)
  
  
  # Create results data.frame
  results <- data.table(
    Model = names(models),
    AIC = aic_values,
    Delta_AIC = delta_aic,
    BIC = bic_values,
    Delta_BIC = delta_bic,
    row.names = NULL
  )
  
  # Sort by AIC
  results <- results[order(results$AIC), ]
  
  # Add model descriptions
  model_descriptions <- c(
    full = "Full model (latitude for mu, sigma, and nu)",
    zero_only = "Zero-inflation only (latitude for nu, intercept for mu)",
    magnitude_only = "Magnitude only (latitude for mu, intercept for nu)",
    intercept_only = "Partial null (latitude for sigma, intercepts for mu and nu)",
    complete_null = "Complete null (intercepts only)"
  )
  
  results$Description <- model_descriptions[results$Model]
  
  cli::cli_alert_success("Model fitting finished successfully")
  
  # Return both the results table and the list of models
  return(list(
    summary = results,
    models = models
  ))
}

print_gamlss_summary <- function(models, results) {
  # Print results with cli
  cli_h2("Model Comparison Results")
  cli_text("")
  
  # Create a rule for the table header
  cli_rule(left = col_br_blue("Model Comparison Table"))
  
  # Function to get colored delta value
  get_colored_delta <- function(delta) {
    if (delta < 2) {
      return(cli::col_green(sprintf("%.2f", delta)))
    } else if (delta < 4) {
      return(cli::col_yellow(sprintf("%.2f", delta)))
    } else {
      return(cli::col_red(sprintf("%.2f", delta)))
    }
  }
  
  # Print each model result
  for (i in 1:nrow(results)) {
    model_name <- results$Model[i]
    aic_value <- results$AIC[i]
    bic_value <- results$BIC[i]
    delta_aic <- results$Delta_AIC[i]
    delta_bic <- results$Delta_BIC[i]
    description <- results$Description[i]
    
    # Color coding based on delta AIC
    if (delta_aic < 2 && delta_bic < 2) {
      status_symbol <- col_green(cli::symbol$tick)
    } else if (delta_aic < 4 || delta_bic < 4) {
      status_symbol <- make_ansi_style("#FFA500")(cli::symbol$warning)
    } else {
      status_symbol <- col_red(cli::symbol$cross)
    }
    
    # Get model formulas
    model <- models[[model_name]]
    mu_formula <- deparse(model$mu.formula)
    sigma_formula <- deparse(model$sigma.formula)
    nu_formula <- deparse(model$nu.formula)
    
    # Get fit statistics
    df_fit <- model$df.fit
    df_residual <- model$df.residual
    global_dev <- model$G.deviance
    
    # Print model information
    cli_alert(paste(status_symbol,
                    "{.field {model_name}} ",
                    "{.emph {description}}", collapse = " "))
    
    # Print formulas with proper formatting
    cli_bullets(c(
      " " = cli::col_cyan(paste0("μ: ", mu_formula)),
      " " = cli::col_cyan(paste0("σ: ", sigma_formula)),
      " " = cli::col_cyan(paste0("ν: ", nu_formula))
    ))
    
    
    # Print detailed fit statistics
    cli_bullets(c(
      " " = paste0("Fit Statistics:"),
      " " = paste0("Parameters (df): ", df_fit),
      " " = paste0("Residual df: ", df_residual),
      " " = paste0("Global Deviance: ", round(global_dev, 2)),
      " " = paste0("AIC: ", round(aic_value, 2), "  (ΔAIC: ", get_colored_delta(delta_aic), ")"),
      " " = paste0("BIC/SBC: ", round(bic_value, 2), "  (ΔBIC: ", get_colored_delta(delta_bic), ")")
    ))
    
    cli_text("")
    cli_text("")
  }
  
  # Print interpretation
  cli_h2("Interpretation")
  best_aic_model <- results$Model[which.min(results$AIC)]
  best_bic_model <- results$Model[which.min(results$BIC)]
  
  if (best_aic_model == best_bic_model) {
    cli_alert_success("Both AIC and BIC support the {.strong {best_aic_model}} model")
  } else {
    cli_alert_warning("AIC supports {.strong {best_aic_model}} while BIC supports {.strong {best_bic_model}}")
    cli_text("")
    cli_text("Note: BIC tends to favor simpler models, especially with larger sample sizes")
  }
}

compare_gam_models <- function(data, predictor = "centroid_latitude", response = "overlapRegion") {
  # Create model formulas with smooth terms
  full_formula <- as.formula(paste(response, "~ s(", predictor, ")"))
  null_formula <- as.formula(paste(response, "~ 1"))
  
  # Start CLI output
  cli_h1("GAM Model Comparison")
  cli_alert_info("Fitting models...")
  
  # List to store models
  models <- list()
  
  # Progress bar
  cli_progress_bar("Fitting models", total = 2)
  
  # 1. Full model with smooth term
  cli_progress_update()
  models$full <- gam(
    formula = full_formula,
    family = betar(link = "logit"),
    data = data,
    method = "REML"  # Restricted Maximum Likelihood for smooth parameter estimation
  )
  
  # 2. Null model
  cli_progress_update()
  models$null <- gam(
    formula = null_formula,
    family = betar(link = "logit"),
    data = data,
    method = "REML"
  )
  
  cli_progress_done()
  
  # Calculate AIC and BIC
  aic_values <- sapply(models, AIC)
  bic_values <- sapply(models, BIC)
  delta_aic <- aic_values - min(aic_values)
  delta_bic <- bic_values - min(bic_values)
  
  # Create results data.table
  results <- data.table(
    Model = names(models),
    AIC = aic_values,
    Delta_AIC = delta_aic,
    BIC = bic_values,
    Delta_BIC = delta_bic,
    row.names = NULL
  )
  
  # Sort by AIC
  results <- results[order(results$AIC), ]
  
  # Add model descriptions
  model_descriptions <- c(
    full = "Full model (smooth term for latitude)",
    null = "Null model (intercept only)"
  )
  
  results$Description <- model_descriptions[results$Model]
  
  cli::cli_alert_success("Model fitting finished successfully")
  
  return(list(
    summary = results,
    models = models
  ))
}
