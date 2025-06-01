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
    vebprint(unique(expected_names), verbose, text = "Expected:")

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

  vebcat("Stats successfully initialised.", color = "funSuccess")

  return(list(
    included = included_sp,
    excluded = excluded_sp
  ))
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
  init <- function(spec.filename, region, extra, verbose) {
    sp_name <- clean_spec_filename(dirname(spec.filename[1]))
    vebcat("Species:", sp_name)

    sp_rast <- load_sp_rast(spec.filename[1])
    vebprint(sp_rast, text = "Raster:")

    summed_dt <- extract_raster_to_dt(sp_rast, value = "cellRichness", cells = TRUE)
    summed_dt <- unique(summed_dt, by = "cell")

    return(summed_dt)
  }

  execute <- function(spec.filename, region, extra, verbose) {
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
        fifelse(
          is.na(richnessSum) & is.na(i.cellAbundance), NA_real_,
          fifelse(
            is.na(richnessSum), i.cellAbundance,
            fifelse(
              is.na(i.cellAbundance), richnessSum,
              richnessSum + i.cellAbundance
            )
          )
        )]

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

      dt_summed[dt_sp,
        on = "cell", cellRichness :=
          fifelse(
            is.na(cellRichness) & is.na(i.richnessSum), NA_real_, # NA + NA = NA
            fifelse(
              is.na(cellRichness), i.richnessSum, # NA + value = value
              fifelse(
                is.na(i.richnessSum), cellRichness, # value + NA = value
                cellRichness + i.richnessSum
              )
            )
          ) # value1 + value2 = sum of values
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
  vebcat("Setting up world Map.", veb = verbose, color = "funInit")

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

  world_map <- suppressWarnings(check_crs(world_map, projection, verbose))

  vebcat("World Map setup completed successfully", veb = verbose, color = "funSuccess")

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
    # sp_regions_dt <- merge(sp_dt, region_dt, by = "ID")
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
      country = first(country), # Assuming country is the same for each subRegionName
      species = first(species), # Assuming species is the same for all rows
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
        source("./R/utils/utils.R")
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
        tryCatch(
          {
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
          }
        )
      })
    }, error = function(e) {
      vebcat("An error occurred in the parallel process ~ stopping cluster and closing connections.", color = "fatalError")
      catn("Error messages written to node logs in:\n", colcat(nodes_dir, color = "output"))
      stop(e)
    }, finally = {
      catn("Cleaning up cluster connections")
      if (exists("cl")) {
        tryCatch(
          {
            stopCluster(cl)
          },
          error = function(e) {
            vebcat("Error stopping cluster:", e$message, color = "fatalError")
          }
        )
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

calculate_taxon_richness <- function(dt, taxon, region = TRUE, country = NULL, verbose = FALSE) {
  catn("Calculating richness for", highcat(taxon))

  by_if_region <- if (region) c(taxon, "subRegionName") else c(taxon)
  by_taxon_or_region <- if (region) "subRegionName" else NULL

  vebprint(by_if_region, text = "by_if_region:")
  vebprint(by_taxon_or_region, text = "by_taxon_or_region:")
  vebprint(region, text = "region:")

  vebprint(ifelse(region, "paoo", "cleanName"), text = "ifelse:")

  dt_copy <- copy(dt)
  # Calculate richness per taxon and subregion (number of unique species in each combination)
  dt_copy[, taxonRichness := ifelse(region,
    length(unique(cleanName[PAoO > 0])),
    length(unique(cleanName))
  ),
  by = by_if_region
  ]
  # How many unique species are in each order

  vebprint(dt_copy, verbose, text = "Taxon Richness before sum:")

  # Calculate total richness per subregion (across all taxa)
  dt_copy[, totalRichness := ifelse(region,
    length(unique(cleanName[PAoO > 0])),
    .N
  ),
  by = by_taxon_or_region
  ]

  # Calculate relative richness
  dt_copy[, relativeRichness := taxonRichness / totalRichness, by = by_if_region]

  # Remove duplicates to get one row per taxon-subregion combination
  if (is.null(country)) {
    dt_copy <- unique(dt_copy, by = by_if_region)
  } else {
    by_if_region <- c(by_if_region, country)
    dt_copy <- unique(dt_copy, by = by_if_region)
  }

  vebprint(dt_copy, verbose, text = "Taxon Richness after sum:")

  return(dt_copy)
}

get_taxon_richness <- function(paoo.file, stats, taxon, country = FALSE, verbose = FALSE) {
  vebcat("Getting region richness", color = "funInit")

  if (taxon == "species") taxon <- "cleanName"

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


get_unknown_composition <- function(unknown.path, stats, taxon = "order", out.dir, verbose = FALSE) {
  out_file <- file.path(out.dir, "unknown-composition.csv")

  if (file.exists(out_file)) {
    richness_dt <- fread(out_file)
  } else {
    unknown_dt <- fread(unknown.path, sep = "\t")

    unknown_dt[, c("genus", "specificEpithet", "cleanName", "other", "fullName", "extra", "structure") := {
      res <- lapply(seq_along(scientificName), function(i) {
        cat("\rProcessing rows for scientificName:", i, "of", .N)
        flush.console()
        clean_spec_name(scientificName[i], config$species$standard_symbols, config$species$standard_infraEpithets, verbose)
      })
      catn()
      list(
        vapply(res, function(x) x$genus, character(1)),
        vapply(res, function(x) x$specificEpithet, character(1)),
        vapply(res, function(x) x$cleanName, character(1)),
        vapply(res, function(x) x$other, character(1)),
        vapply(res, function(x) x$fullName, character(1)),
        vapply(res, function(x) x$extra, character(1)),
        vapply(res, function(x) x$structure, character(1))
      )
    }]

    setnames(unknown_dt, "scientificName", "verbatim_name")

    # Get the taxonomic hierarchy
    unknown_taxa <- as.data.table(get_spec_taxons(unknown_dt$verbatim_name))

    # find the taxon group
    richness_dt <- get_order_group(
      unknown_taxa,
      spec.col = "species"
    )

    # Merge
    richness_dt <- richness_dt[unknown_dt, on = "verbatim_name", nomatch = NULL]

    # Remove columns with names containing "i."
    richness_dt <- richness_dt[, !grepl("^i\\.", names(richness_dt)), with = FALSE]

    # Calculate richness
    richness_dt <- calculate_taxon_richness(
      dt = richness_dt,
      taxon = taxon,
      region = FALSE,
      verbose = verbose
    )

    fwrite(richness_dt, out_file, bom = TRUE)
  }


  return(richness_dt)
}

# Function to calculate and export differences between study data and GloNAF
calculate_order_differences <- function(dt, dt_comparison, vis.x, vis.y, vis.fill, out.dir) {
  out_file <- file.path(out.dir, "order_differences.csv")

  if (file.exists(out_file)) {
    result <- fread(out_file)
  } else {
    # Create copy of data
    dt_copy <- copy(dt)
    comparison_copy <- copy(dt_comparison)

    # Get all unique orders
    all_orders <- unique(c(dt_copy[[vis.fill]], comparison_copy[[vis.fill]]))

    # Calculate total richness for each region and for GloNAF
    region_totals <- dt_copy[, .(total_richness = sum(get(vis.y))), by = vis.x]
    glonaf_total <- sum(comparison_copy[[vis.y]])

    # Calculate GloNAF relative richness for each order
    glonaf_relative <- comparison_copy[, .(glonaf_relative = sum(get(vis.y)) / glonaf_total), by = vis.fill]

    # Calculate relative richness for each region and order
    region_order_richness <- dt_copy[, .(study_relative = sum(get(vis.y)) / region_totals[get(vis.x) == .BY[[vis.x]], total_richness]),
      by = c(vis.x, vis.fill)
    ]

    # Merge with GloNAF relative richness
    result <- region_order_richness[glonaf_relative, on = vis.fill, all = TRUE]

    # Replace NA with 0
    result[is.na(study_relative), study_relative := 0]
    result[is.na(glonaf_relative), glonaf_relative := 0]

    # Remove NA orders
    result <- result[!is.na(get(vis.fill))]

    # Calculate difference
    result[, difference := study_relative - glonaf_relative]

    # Add absolute difference column for sorting
    result[, abs_difference := abs(difference)]

    # Sort by absolute difference
    setorder(result, -abs_difference)

    # Save to CSV
    fwrite(result, out_file)
  }

  return(result)
}

analyze_composition_significance <- function(study_dt, comparison_dt, taxon_col = "order", region_col = "subRegionName", richness_col = "relativeRichness", alpha = 0.05, save.dir = NULL, save.device = "jpeg", plot.show = FALSE, verbose = FALSE, min_sample_size = 30) {
  vebcat("Testing taxonomic composition significance", color = "funInit")
  
  # Get GloNAF baseline proportions
  glonaf_props <- comparison_dt[, .(
    glonaf_count = sum(get(richness_col))
  ), by = taxon_col]
  
  glonaf_total <- sum(glonaf_props$glonaf_count)
  glonaf_props[, expected_prop := glonaf_count / glonaf_total]
  
  vebprint(glonaf_props, verbose, "GloNAF baseline proportions:")
  
  # Get unique regions
  regions <- unique(study_dt[[region_col]])
  
  # Initialize results table
  results <- data.table(
    region = character(),
    chi_square = numeric(),
    df = integer(),
    p_value = numeric(),
    significant = logical(),
    effect_size = numeric(),
    interpretation = character()
  )
  
  # Test each region
  for (region in regions) {
    cat("\rTesting region:", region)
    
    # Get observed data for this region
    region_data <- study_dt[get(region_col) == region]
    
    # Calculate observed proportions (already relative richness)
    region_total <- sum(region_data[[richness_col]], na.rm = TRUE)
    
    # Skip if no data
    if (region_total == 0 || nrow(region_data) == 0) {
      vebcat("No data for", region, color = "warning")
      results <- rbind(results, data.table(
        region = region,
        chi_square = NA,
        df = NA,
        p_value = NA,
        significant = NA,
        effect_size = NA,
        interpretation = "No data"
      ), fill = TRUE)
      next
    }
    
    # Normalize to proportions if not already
    region_data[, observed_prop := get(richness_col) / region_total]
    
    # Get all taxonomic orders present in either dataset
    all_orders <- unique(c(region_data[[taxon_col]], glonaf_props[[taxon_col]]))
    
    # Create complete data frame with zeros for missing orders
    complete_data <- data.table(order = all_orders)
    setnames(complete_data, "order", taxon_col)
    
    # Join with observed data using data.table syntax
    complete_data <- region_data[, c(taxon_col, "observed_prop"), with = FALSE][complete_data, on = taxon_col]
    complete_data[is.na(observed_prop), observed_prop := 0]
    
    # Join with expected proportions using data.table syntax
    complete_data <- glonaf_props[, c(taxon_col, "expected_prop"), with = FALSE][complete_data, on = taxon_col]
    complete_data[is.na(expected_prop), expected_prop := 0]
    
    # Remove orders with zero expected proportions (can't test these)
    test_data <- complete_data[expected_prop > 0]
    
    # Convert proportions to counts for chi-square test using artificial sample size
    test_data[, observed_count := round(observed_prop * min_sample_size)]
    test_data[, expected_count := round(expected_prop * min_sample_size)]
    
    # Ensure we have at least some observations
    total_obs <- sum(test_data$observed_count)
    if (total_obs == 0) {
      # If all zeros, assign at least 1 count to maintain proportions
      max_prop_idx <- which.max(test_data$observed_prop)
      if (length(max_prop_idx) > 0) {
        test_data[max_prop_idx, observed_count := 1]
        total_obs <- 1
      }
    }
    
    # Check if we have enough data for chi-square test
    if (nrow(test_data) < 2) {
      vebcat("Insufficient taxonomic diversity for", region, color = "warning")
      
      results <- rbind(results, data.table(
        region = region,
        chi_square = NA,
        df = NA,
        p_value = NA,
        significant = NA,
        effect_size = NA,
        interpretation = "Insufficient taxonomic diversity"
      ), fill = TRUE)
      next
    }
    
    # Remove categories with zero expected counts
    test_data <- test_data[expected_count > 0]
    
    if (nrow(test_data) < 2) {
      results <- rbind(results, data.table(
        region = region,
        chi_square = NA,
        df = NA,
        p_value = NA,
        significant = NA,
        effect_size = NA,
        interpretation = "No matching taxonomic orders"
      ), fill = TRUE)
      next
    }
    
    # Perform chi-square goodness-of-fit test
    tryCatch(
      {
        # Use the proportions directly, not counts, for the test
        chi_test <- chisq.test(
          x = test_data$observed_count,
          p = test_data$expected_prop / sum(test_data$expected_prop), # Renormalize
          simulate.p.value = TRUE, # Always use Monte Carlo for reliability
          B = 2000
        )
        
        # Calculate effect size (Cramér's V)
        n <- sum(test_data$observed_count)
        v <- sqrt(chi_test$statistic / (n * (length(test_data$expected_prop) - 1)))
        
        # Interpret effect size using data.table syntax
        interpretation <- fifelse(is.na(v), "Cannot calculate",
                                  fifelse(v < 0.1, "Negligible difference",
                                          fifelse(v < 0.3, "Small difference", 
                                                  fifelse(v < 0.5, "Medium difference", 
                                                          "Large difference"))))
        
        # Add to results
        results <- rbind(results, data.table(
          region = region,
          chi_square = as.numeric(chi_test$statistic),
          df = chi_test$parameter,
          p_value = chi_test$p.value,
          significant = chi_test$p.value < alpha,
          effect_size = as.numeric(v),
          interpretation = interpretation
        ), fill = TRUE)
      },
      error = function(e) {
        vebcat("Error testing", region, ":", e$message, color = "warning")
        
        results <- rbind(results, data.table(
          region = region,
          chi_square = NA,
          df = NA,
          p_value = NA,
          significant = NA,
          effect_size = NA,
          interpretation = paste("Test failed:", e$message)
        ), fill = TRUE)
      }
    )
  }
  
  catn()
  
  # Summary statistics
  valid_tests <- results[!is.na(p_value)]
  n_significant <- sum(valid_tests$significant, na.rm = TRUE)
  n_total <- nrow(valid_tests)
  
  # Apply multiple testing correction
  if (nrow(valid_tests) > 1) {
    results[!is.na(p_value), p_adjusted := p.adjust(p_value, method = "fdr")]
    results[!is.na(p_adjusted), significant_adjusted := p_adjusted < alpha]
    n_significant_adj <- sum(results[!is.na(p_adjusted), p_adjusted < alpha], na.rm = TRUE)
  } else {
    results[, p_adjusted := p_value]
    results[, significant_adjusted := significant]
    n_significant_adj <- n_significant
  }
  
  # Print summary
  cli::cli_h2("Taxonomic Composition Significance Test Results")
  cli::cli_bullets(c(
    "*" = "Total regions tested: {.field {n_total}}",
    "*" = "Regions with significant differences (uncorrected): {.field {n_significant}} ({round(100*n_significant/n_total, 1)}%)",
    "*" = "Regions with significant differences (FDR corrected): {.field {n_significant_adj}} ({round(100*n_significant_adj/n_total, 1)}%)",
    "*" = "Null hypothesis: Taxonomic composition matches GloNAF baseline",
    "*" = "Alternative: Taxonomic composition differs from GloNAF baseline"
  ))
  
  # Sort by effect size for easy interpretation
  results <- results[order(-effect_size)]
  
  # Create visualizations (rest of visualization code remains the same)
  if (!is.null(save.dir)) {
    create_dir_if(save.dir)
    
    vebcat("Creating significance test visualizations", color = "funInit")
    
    # 1. P-values and effect sizes plot
    if (nrow(valid_tests) > 0) {
      p1 <- ggplot(results, aes(x = reorder(region, -effect_size))) +
        geom_col(aes(y = effect_size, fill = significant_adjusted), alpha = 0.7) +
        geom_hline(yintercept = c(0.1, 0.3, 0.5), linetype = "dashed", color = "gray50") +
        scale_fill_manual(
          values = c("TRUE" = "#d73027", "FALSE" = "#4575b4"),
          labels = c("TRUE" = "Significant", "FALSE" = "Not Significant"),
          name = "FDR-corrected"
        ) +
        labs(
          title = "Effect Sizes of Taxonomic Composition Differences",
          subtitle = "Cramér's V comparing each region to GloNAF baseline",
          x = "Arctic Subregion",
          y = "Effect Size (Cramér's V)",
          caption = "Dashed lines: 0.1=small, 0.3=medium, 0.5=large effect"
        ) +
        theme_minimal() +
        ggplot.theme() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          legend.position = "bottom"
        )
      
      save_ggplot(
        save.plot = p1,
        save.name = "composition-significance-overview",
        save.width = 3000,
        save.height = 2400,
        save.dir = save.dir,
        save.device = save.device,
        plot.show = plot.show,
        verbose = verbose
      )
    }
    
    # Rest of visualization code continues...
    vebcat("Significance test visualizations completed", color = "funSuccess")
  }
  
  vebcat("Composition significance testing completed", color = "funSuccess")
  
  return(results)
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
    (bbox[2] - bbox[1]) * 0.1, # 10% padding
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
  if (aspect_ratio > 1) { # wider than tall
    base_width <- plot_x_range * width
    base_height <- plot_y_range * height
    actual_width <- base_width
    actual_height <- base_width / aspect_ratio
  } else { # taller than wide
    base_width <- plot_x_range * (width + 0.1)
    base_height <- plot_y_range * (height + 0.1)

    actual_width <- base_height * aspect_ratio
    actual_height <- base_height
  }

  # Position consistently from top edge
  inset_ymax <- unname(plot.ext[4]) - (top * actual_height)
  inset_ymin <- inset_ymax - base_height # Use base_height for consistent top spacing

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

    cleaned_data[, `:=`(
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
    nu.formula = ~1,
    family = BEZI,
    data = data,
    trace = FALSE
  )

  # 4. Intercept only model (intercepts for both nu and mu)
  cli_progress_update()
  models$intercept_only <- gamlss(
    formula = null_formula,
    sigma.formula = as.formula(paste("~", predictor)),
    nu.formula = ~1,
    family = BEZI,
    data = data,
    trace = FALSE
  )

  # 5. Complete null model (intercepts for mu, sigma, and nu)
  cli_progress_update()
  models$complete_null <- gamlss(
    formula = null_formula,
    sigma.formula = ~1,
    nu.formula = ~1,
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
      "{.emph {description}}",
      collapse = " "
    ))

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

compare_gam_models <- function(data, predictor = "centroidLatitude", response = "overlapRegion") {
  if (any(data[[response]] < 0 | data[[response]] > 1)) {
    cli::cli_alert_danger("Response variable must be between 0 and 1 exclusively for beta regression")
    stop("Invalid response variable range")
  }

  # Create model formulas with smooth terms
  full_formula <- as.formula(paste(response, "~ s(", predictor, ")"))
  null_formula <- as.formula(paste(response, "~ 1"))

  # Start CLI output
  cli_h1("GAM Model Comparison")
  cli_alert_info("Fitting models...")

  # List to store models
  models <- list()

  # Progress bar
  cli_progress_bar("Fitting models", total = 2, type = "download")

  # 1. Full model with smooth term
  cli_progress_update()
  models$full <- gam(
    formula = full_formula,
    family = betar(link = "logit"),
    data = data,
    method = "REML" # Restricted Maximum Likelihood for smooth parameter estimation
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

  # Get scores
  reml_scores <- sapply(models, function(x) summary(x)$sp.criterion)
  delta_reml <- reml_scores - min(reml_scores)
  precision_scores <- sapply(models, function(x) summary(x)$family$getTheta())
  r_sq_scores <- sapply(models, function(x) summary(x)$r.sq)
  deviance_scores <- sapply(models, function(x) summary(x)$dev.expl)
  edf_scores <- c(sum(summary(models$full)$s.table[, 1]), 0)
  edf_total_scores <- c(sum(summary(models$full)$s.table[, 2]), 0)

  results <- data.table::data.table(
    model = names(models),
    reml = reml_scores,
    delta_reml = delta_reml,
    precision = precision_scores,
    edf = edf_scores,
    edf.total = edf_total_scores,
    r_sq_adj = r_sq_scores,
    deviance_explained = deviance_scores,
    description = c(
      full = "Full model (smooth term for latitude)",
      null = "Null model (intercept only)"
    )
  )
  # Sort by REML score (lower is better)
  results <- results[order(results$reml), ]

  cli::cli_alert_success("Model fitting finished successfully")

  result <- list(
    summary = results,
    models = models
  )

  print_gam_summary(models = result$models, results = result$summary)

  return(result)
}

print_gam_summary <- function(models, results) {
  # Helper function for consistent number formatting
  format_num <- function(x, digits = 3) {
    format(round(as.numeric(x), digits), big.mark = ",", trim = TRUE)
  }

  # Color schemes for different elements
  heading_color <- "cyan"
  subheading_color <- "magenta"
  stat_color <- "green"
  warning_color <- "yellow"

  # ---- Model Specifications Section ----
  cli::cli_h1(cli::col_cyan("GAM Model Analysis"))

  for (i in 1:length(models)) {
    model_summary <- summary(models[[i]])
    model_name <- names(models)[i]

    # Model header with colored name
    cli::cli_h2(cli::col_magenta(paste("Model:", model_name)))

    # Core model statistics with colored values
    cli::cli_bullets(c(
      "►" = sprintf(
        "Family: %s {.val %s} (Precision: {.val %.3f})",
        models[[i]]$family$family,
        models[[i]]$family$link,
        models[[i]]$family$getTheta()
      ),
      "►" = sprintf("R-sq.(adj): {.val %.4f}", model_summary$r.sq),
      "►" = sprintf("Deviance explained: {.val %.2f%%}", model_summary$dev.expl * 100),
      "►" = sprintf("REML score: {.val %s}", format_num(model_summary$sp.criterion)),
      "►" = sprintf("Scale est.: {.val %.3f}", model_summary$scale),
      "►" = sprintf("Sample size: {.val %d}", model_summary$n)
    ))

    # ---- Coefficients Tables ----
    # Parametric coefficients
    cli::cli_h3(cli::col_cyan("Parametric Coefficients"))
    if (is.data.table(model_summary$p.table)) {
      print_formatted_table(model_summary$p.table, "parametric")
    } else {
      print(model_summary$p.table)
    }

    # Smooth terms
    cli::cli_h3(cli::col_cyan("Smooth Terms"))
    if (is.data.table(model_summary$s.table)) {
      print_formatted_table(model_summary$s.table, "smooth")
    } else {
      print(model_summary$s.table)
    }

    cli::cli_rule()
  }

  # ---- Model Comparison Section ----
  cli::cli_h2(cli::col_magenta("Model Comparison"))

  # Format results columns
  results[, `:=`(
    REML = format_num(reml, 1),
    Delta_REML = format_num(delta_reml, 1),
    Precision = format_num(precision, 3),
    EDF = format_num(edf, 2),
    EDF.Total = format_num(edf.total, 2),
    R_sq_adj = format_num(r_sq_adj, 4),
    Deviance_Explained = paste0(format_num(deviance_explained * 100, 2), "%")
  )]

  # Print comparison statistics with enhanced formatting
  cli::cli_h3(cli::col_cyan("Model Performance Metrics"))
  cli::cli_bullets(c(
    "►" = cli::col_green(sprintf(
      "REML Difference: {.val %s} (lower is better)",
      format_num(abs(as.numeric(gsub(",", "", results$REML[1])) -
        as.numeric(gsub(",", "", results$REML[2]))), 1)
    )),
    "►" = cli::col_green(sprintf(
      "Model Complexity (EDF): Full = {.val %s}, Null = {.val %s}",
      results$EDF[1], results$EDF[2]
    )),
    "►" = cli::col_green(sprintf(
      "Precision Parameters: Full = {.val %s}, Null = {.val %s}",
      results$Precision[1], results$Precision[2]
    )),
    "►" = cli::col_green(sprintf(
      "R-sq.(adj): Full = {.val %s}, Null = {.val %s}",
      results$R_sq_adj[1],
      ifelse(length(results$R_sq_adj) > 1, results$R_sq_adj[2], "N/A")
    )),
    "►" = cli::col_green(sprintf(
      "Deviance Explained: Full = {.val %s}, Null = {.val %s}",
      results$Deviance_Explained[1],
      ifelse(length(results$Deviance_Explained) > 1, results$Deviance_Explained[2], "N/A")
    ))
  ))

  # ---- Interpretation Guide ----
  cli::cli_h3(cli::col_magenta("Interpretation Guide"))
  cli::cli_bullets(c(
    "!" = cli::col_yellow("REML scores: Lower values indicate better model fit"),
    "!" = cli::col_yellow("EDF: Higher values show increased model complexity"),
    "!" = cli::col_yellow("Precision: Higher values indicate less data dispersion"),
    "!" = cli::col_yellow("R-sq.(adj): Proportion of variance explained (higher is better)"),
    "!" = cli::col_yellow("Deviance explained: Percentage improvement over null model")
  ))

  invisible(results)
}

# Enhanced table printing helper function
print_formatted_table <- function(table_data, table_type) {
  col_names <- colnames(table_data)

  # Add significance stars
  if (any(c("Pr(>|z|)", "p-value") %in% col_names)) {
    p_col <- ifelse("Pr(>|z|)" %in% col_names, "Pr(>|z|)", "p-value")
    table_data$Signif. <- sapply(table_data[[p_col]], function(p) {
      if (p < 0.001) {
        return(cli::col_green("***"))
      }
      if (p < 0.01) {
        return(cli::col_green("**"))
      }
      if (p < 0.05) {
        return(cli::col_green("*"))
      }
      if (p < 0.1) {
        return(cli::col_yellow("."))
      }
      return("")
    })
    col_names <- c(col_names, "Signif.")
  }

  # Print formatted table
  header <- paste(cli::col_cyan(cli::style_bold(col_names)), collapse = "    ")
  cli::cli_text(header)
  cli::cli_rule(col = "cyan")

  for (i in 1:nrow(table_data)) {
    row_data <- format(table_data[i, ], digits = 4, trim = TRUE)
    row_text <- paste(row_data, collapse = "    ")
    cli::cli_text(row_text)
  }

  cli::cli_rule(col = "cyan")
}


calculate_species_centroids <- function(dt, out.dir, verbose = FALSE) {
  # Input dt should have columns: species, level3Long, level3Lat
  out_file <- file.path(out.dir, "spec-centroids.csv")

  if (file.exists(out_file)) {
    species_centroids <- fread(out_file)
  } else {
    dt[, inRegion := NULL]

    spec_dt <- rbindlist(sp_stats)

    spec_dt[, `:=`(
      decimalLongitude = level3Long,
      decimalLatitude = level3Lat,
      level3Long = NULL,
      level3Lat = NULL,
      level3Code = NULL,
      level2Code = NULL,
      level1Code = NULL
    )]

    spec_centroids <- data.table()

    unique_species <- unique(spec_dt$cleanName)

    catn("Iterations:", highcat(length(unique_species)))

    for (i in cli_progress_along(1:length(unique_species), "downloading")) {
      spec <- spec_dt[cleanName == unique_species[i]]

      spec_out <- find_wgsrpd_region(
        spec.dt = spec,
        projection = "longlat",
        longitude = "decimalLongitude",
        latitude = "decimalLatitude",
        wgsrpd.dir = "./resources/region/wgsrpd",
        wgsrpdlvl = "3",
        wgsrpdlvl.name = TRUE,
        unique = TRUE,
        suppress = TRUE,
        verbose = verbose
      )

      spec_out <- spec[spec_out, on = "level3Name"]

      spec_centroids <- rbind(spec_centroids, spec_out)
    }

    fwrite(spec_centroids, out_file, bom = TRUE)
  }

  return(species_centroids)
}

prepare_poly_spec <- function(dt, response, predictor, transform = NULL, by.region = FALSE) {
  vebcat("Preparing polynomial species", color = "funInit")
  # Keep one row per species per level3 botanical country
  dt <- dt[, .SD, .SDcols = c(response, predictor), by = .(cleanName, level3Name)]

  if (by.region) {
    dt <- dt[, .(
      res = mean(get(response), na.rm = TRUE),
      pred = first(get(predictor))
    ), by = level3Name]

    # Rename the columns to match the input names
    setnames(dt, c("res", "pred"), c(response, predictor))
  }

  dt <- dt[!is.na(dt[[response]]) & !is.na(dt[[predictor]])]

  dt[get(response) > 1, (response) := 1]
  dt[get(response) < 0, (response) := 0]

  if (!is.null(transform) && transform == "logit") {
    logit <- function(p) {
      p <- ifelse(p == 0, 1e-6, ifelse(p == 1, 1 - 1e-6, p))
      log(p / (1 - p))
    }

    # Apply the logit transformation to the response variable
    dt[[response]] <- logit(dt[[response]])
  } else {
    if (any(dt[[response]] < 0 | dt[[response]] > 1)) {
      stop("Response variable must be between 0 and 1")
    }
  }

  setorderv(dt, predictor)

  return(dt)
}

fit_polynomial_model <- function(dt, response = "overlapRegion", predictor = "centroidLat") {
  vebcat("Fitting polynomial model", color = "funInit")

  # Create B-spline basis with knot at equator (0 degrees)
  # Reduce degrees of freedom for smoother curves by setting df parameter
  bs_matrix <- bs(dt[[predictor]],
    knots = 0,
    # df = 10,      # Control smoothness - fewer df = smoother curves
    degree = 3
  ) # Cubic spline

  # Define quantiles to model
  taus <- c(0.1, 0.25, 0.5, 0.75, 0.9)

  vebcat("Fitting models for each quantile")
  # Fit models for each quantile using logistic transformation
  models <- lapply(taus, function(tau) {
    rq(as.formula(paste(response, "~ bs_matrix")), data = dt, tau = tau, method = "br")
  })

  # Create prediction grid
  pred_grid <- data.table(
    lat = seq(min(dt[[predictor]]), max(dt[[predictor]]), length.out = nrow(bs_matrix))
  )

  vebcat("Acquiring predictions for each quantile")
  # Get predictions for each quantile
  pred_bs <- bs(pred_grid$lat,
    knots = 0,
    degree = 3,
    Boundary.knots = range(dt[[predictor]])
  )

  pred_dt <- as.data.table(pred_bs)

  # Get predictions for each quantile
  predictions <- sapply(models, function(model) {
    predict(model, newdata = pred_dt)
  })

  # Combine predictions with grid
  pred_dt <- data.table(
    latitude = pred_grid$lat,
    q10 = predictions[, 1],
    q25 = predictions[, 2],
    q50 = predictions[, 3],
    q75 = predictions[, 4],
    q90 = predictions[, 5]
  )[order(latitude)]

  return(list(
    models = models,
    predictions = pred_dt,
    taus = taus,
    bs_matrix = bs_matrix
  ))
}

summarize_polynomial_model <- function(model_results, transform = NULL) {
  # Extract coefficients for each quantile
  cli::cli_h1("B-spline Quantile Regression Model Summary")

  # Get coefficients for each quantile
  cli::cli_h2("Model Coefficients")
  coef_summary <- lapply(seq_along(model_results$models), function(i) {
    tau <- model_results$taus[i]
    model <- model_results$models[[i]]
    coefs <- coef(model)
    if (!is.null(transform)) odds_ratios <- exp(coefs)

    # Print coefficient summary for each quantile
    cli::cli_alert_info("Quantile {.val {tau * 100}}%")
    coef_table <- data.table(
      term = c(
        "intercept",
        "cubic_B-spline_basis1", # or spline_basis1
        "cubic_B-spline_basis2",
        "cubic_B-spline_basis3",
        "cubic_B-spline_basis4"
      ),
      estimate = coefs,
      odd_ratios = if (!is.null(transform)) odds_ratios else NULL
    )

    # Format coefficient display
    cli::cli_text("")
    print(coef_table)
    cli::cli_text("")

    data.table(
      quantile = paste0("q", tau * 100),
      term = names(coefs),
      estimate = coefs,
      odds_ratio = if (!is.null(transform)) odds_ratios else NULL
    )
  })

  coef_summary <- rbindlist(coef_summary)

  # Calculate and display R1 statistics
  cli::cli_h2("Model Fit Statistics")
  r1_stats <- lapply(seq_along(model_results$models), function(i) {
    tau <- model_results$taus[i]
    model <- model_results$models[[i]]

    # Get the response values
    y <- model$y

    # Calculate proportion of zeros
    prop_zeros <- mean(y < 0.001)

    # Calculate pseudo R² for quantile regression
    residuals <- model$residuals
    y_tau <- quantile(y, tau)
    null_deviance <- sum(abs(y - y_tau))
    model_deviance <- sum(abs(residuals))

    # Calculate pseudo R²
    r1 <- 1 - model_deviance / null_deviance

    cli::cli_alert_info(sprintf(
      "Quantile %.0f%%: %.1f%% of observations are near zero",
      tau * 100, prop_zeros * 100
    ))

    if (is.finite(r1)) {
      cli::cli_alert_success("Pseudo R² at {.val {tau * 100}}% quantile: {.val {round(r1, 3)}}")
    }

    data.table(
      quantile = paste0("q", tau * 100),
      R1 = if (is.finite(r1)) r1 else NA,
      prop_zero = prop_zeros
    )
  })

  r1_stats <- rbindlist(r1_stats)

  # Display model information
  cli::cli_h2("Model Information")
  cli::cli_bullets(c(
    "*" = "Number of observations: {.val {nrow(model_results$bs_matrix)}}",
    "*" = "Degrees of freedom: {.val {ncol(model_results$bs_matrix)}}",
    "*" = "Knot position: {.val 0} degrees (equator)",
    "*" = "Boundary knots: {.val {round(attr(model_results$bs_matrix, 'Boundary.knots')[1], 2)}} to {.val {round(attr(model_results$bs_matrix, 'Boundary.knots')[2], 2)}} degrees"
  ))

  # Display interpretation guide
  cli::cli_h2("Interpretation Guide")
  cli::cli_bullets(c(
    "!" = "R1 values closer to 1 indicate better model fit",
    "!" = "Coefficients show the effect of latitude on overlap at each quantile",
    "!" = "The B-spline allows for different patterns in Northern and Southern hemispheres",
    "!" = "Quantiles show the distribution bounds of overlap values",
    if (!is.null(transform)) "!" <- "Odds ratios represent the multiplicative effect on the odds of climatic overlap for a one-unit change in latitude"
  ))

  return(invisible())
}

load_paoo <- function(group, group.dirs, shape, log.dir, sp_stats, verbose = FALSE) {
  for (group in names(sp_group_dirs)) {
    group_dirs <- sp_group_dirs[[group]]$filename

    paoo <- parallel_spec_handler(
      spec.dirs = group_dirs,
      shape = shape,
      dir = paste0(log_dir_aoo, "/", group),
      hv.project.method = "0.5-inclusion",
      out.order = "-TPAoO",
      fun = get_paoo,
      verbose = verbose
    )

    paoo_files[[group]] <- paoo
  }

  paoo_file <- combine_groups(paoo_files, "-TPAoO")

  paoo_dt <- richness_dt[PAoO > 0]

  return(paoo_dt)
}
