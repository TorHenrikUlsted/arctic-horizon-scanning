calc_abs_by_col <- function(data, column = "decimalLatitude", method = "median", verbose = FALSE) {
  if (is.character(data)) {
    if (grepl("\\.txt", data)) {
      vect <- readLines(data)
      dt <- as.data.table(as.numeric(vect))
      setnames(dt, "V1", column)
    } else if (grepl("\\.csv", data)) {
      dt <- fread(data)
    }
  } else if (is.data.table(data)) {
    dt <- copy(data)
  }

  if (method == "median") {
    res <- stats::median(abs(dt[[column]]), na.rm = TRUE)
  } else if (method == "mean") {
    res <- mean(abs(dt[[column]]), na.rm = TRUE)
  }

  return(res)
}


exclude_observations <- function(spec.list, dimensions, method = "median", dt.construct, cols.select, out.file, suppress = FALSE, verbose = FALSE) {
  catn("Checking for species to exclude before analysis.")
  
  if (file.exists(out.file)) {
    final_res <- fread(file.out)
  } else {
    fwrite(dt.construct, out.file, bom = TRUE)
    
    for (i in cli_progress_along(1:length(spec.list), "Downloading")) {
      #cat("\rChecking file", i, "/", length(spec.list))
      spec <- spec.list[[i]]
      spec_name <- basename(gsub(".csv", "", spec))
      spec_dt <- fread(spec, select = cols.select)
      
      nobs <- nrow(spec_dt)
      
      excluded <- FALSE
      
      if (log(nobs) <= length(dimensions)) {
        excluded <- TRUE
      } else {
        next
      }
      
      vebprint(spec_dt, verbose, "Spec:")
      
      spec_condensed <- condense_taxons(spec.dt = spec_dt, suppress = TRUE, verbose)
      
      vebprint(spec_condensed, verbose, "Spec:")
      
      cntry_condensed <- find_wgsrpd_region(
        spec.dt = spec_dt,
        projection = "longlat",
        longitude = "decimalLongitude",
        latitude = "decimalLatitude",
        wgsrpd.dir = "./resources/region/wgsrpd",
        wgsrpdlvl = "3",
        wgsrpdlvl.name = TRUE,
        unique = TRUE,
        suppress = suppress,
        verbose = verbose
      )
      
      vebprint(cntry_condensed, verbose, "Spec:")
      
      # Add columns that are given by the hv sequence
      #lat <- calc_abs_by_col(spec_dt, column = "decimalLatitude")
      
      res <- data.table(
        cleanName = gsub(config$species$file_separator, " ", spec_name),
        iteration = 0,
        observations = nobs,
        dimensions = length(dimensions),
        samplesPerPoint = 0,
        randomPoints = 0,
        excluded = excluded,
        jaccard = 0,
        sorensen = 0,
        fracVolumeSpecies = 1,
        fracVolumeRegion = 1,
        realizedNiche = 0,
        overlapRegion = 0,
        includedOverlap = 0
      )
      
      cbmnd_res <- cbind(res, spec_condensed)
      # Duplicate the rows to match the length of cntry_condensed
      rep_cbmnd_res <- cbmnd_res[rep(seq_len(nrow(cbmnd_res)), nrow(cntry_condensed)), ]
      final_res <- cbind(rep_cbmnd_res, cntry_condensed)
      
      fwrite(final_res, out.file, append = TRUE)
    }
  }
  
  return(final_res)
}

gbif_retry <- function(spec, fun = "name_backbone_checklist", retry.max = 3, timeout = 30, verbose = FALSE) {
  result <- NULL
  if (is.character(fun)) fun <- get(fun)
  for (attempt in 1:retry.max) {
    tryCatch(
      {
        if (attempt == retry.max) {
          closeAllConnections()
          gc()
          timeout <- timeout * 2
        }

        # Check if the function accepts a timeout_ms argument
        if ("timeout_ms" %in% names(formals(fun))) {
          result <- do.call(fun, list(spec, timeout_ms = timeout * 1000))
        } else {
          result <- fun(spec)
        }

        break
      },
      error = function(e) {
        if (attempt == retry.max) {
          vebcat("Failed to fetch data after", retry.max, "attempts. Error:", conditionMessage(e), veb = verbose, color = "warning")
        } else {
          vebcat("Attempt", attempt, "failed. Retrying in", 2^attempt, "seconds.", veb = verbose, color = "nonFatalError")
          Sys.sleep(2^attempt) # Exponential backoff
        }
      }
    )
  }

  if (is.null(result)) {
    vebcat(spec, text = "spec:")
    stop("Failed to fetch data after all retry attempts.")
  }

  return(result)
}

# Find the occurrence points inside the region
occ_region_overlap <- function(spec.occ, region, region.name = "region", spec = "species", long = "decimalLongitude", lat = "decimalLatitude", verbose = FALSE) {
  cols_to_keep <- c(spec, "scientificName", "taxonRank", "speciesKey", long, lat)

  if (is.data.table(spec.occ)) {
    spec_dt <- spec.occ[, ..cols_to_keep]
  } else if (is.character(spec.occ)) {
    catn("Reading data")
    spec_dt <- fread(spec.occ, select = cols_to_keep)
  }

  # Filter out rows with invalid coordinates
  spec_dt <- spec_dt[!is.na(get(long)) & !is.na(get(lat))]

  if (is.character(region)) {
    region <- load_region(region)
  }

  region <- check_crs(region, "longlat")

  # Calculate totalNobs and initialize regionNobs
  spec_dt[, `:=`(
    totalNobs = .N,
    regionNobs = 0L
  ), by = scientificName]

  # Make into points
  catn("Converting to points")
  config_crs <- get_crs_config("longlat")
  points <- terra::vect(spec_dt, geom = c(long, lat), crs = config_crs)

  # Check if overlap with region
  catn("Estimating overlap")
  overlap <- as.data.table(terra::intersect(points, region))

  overlap <- overlap[, .(regionNobs = .N), by = "scientificName"]

  # Merge overlap values with spec data
  spec_dt[overlap, regionNobs := i.regionNobs, on = "scientificName"]

  # Calculate proportional observations
  spec_dt[, propNobs := formatC(regionNobs / totalNobs, format = "e", digits = 3), by = "scientificName"]

  vebcat(highcat(nrow(overlap)), "/", highcat(nrow(spec_dt)), "occurrences found within the region")

  return(spec_dt)
}

loop_occ_overlap <- function(spec.occ.dir, region, region.name = "Arctic", file.out, verbose = FALSE) {
  if (file.exists(file.out)) {
    return(fread(file.out))
  }

  time <- start_timer("Region overlap")

  cols <- c("cleanName", "species", "scientificName", "taxonRank", "speciesKey", "decimalLongitude", "decimalLatitude")
  
  if (file.info(spec.occ.dir[[1]])$isdir) {
    spec_occ <- list.files(spec.occ.dir, full.names = TRUE)
  } else {
    spec_occ <- spec.occ.dir
  }

  # Store unique scientific names for all files
  name_storage <- unique(unlist(lapply(seq_along(spec_occ), function(i) {
    cat("\rStoring species in file:", i, "/", length(spec_occ))
    flush.console()
    unique(fread(spec_occ[i], select = "scientificName"))$scientificName
  })))
  catn()

  # Get name backbone checklist results for all scientificNames
  checklist <- gbif_retry(name_storage, "name_backbone_checklist")
  checklist <- as.data.table(checklist)[, .(species, scientificName, usageKey, status, synonym)]

  print(object.size(name_storage), units = "Mb")
  print(object.size(checklist), units = "Mb")

  rm(name_storage)
  invisible(gc())

  spec_dt_out <- data.table()

  for (i in 1:length(spec_occ)) {
    spec <- spec_occ[i]
    spec_name <- gsub(config$species$file_separator, " ", basename(spec))

    spec_dt <- fread(spec, select = cols)

    cat("\rRunning for species", highcat(i), "/", highcat(length(spec_occ)))
    flush.console()

    # thinned_dt <- clean_spec_occ(
    #   dt = spec_dt,
    #   column = "cleanName",
    #   projection = config_crs <- get_crs_config("longlat"),
    #   resolution = config$projection$raster_scale_m,
    #   seed = config$simulation$seed,
    #   long = "decimalLongitude",
    #   lat = "decimalLatitude",
    #   verbose
    # )

    invisible(capture.output(
      {
        spec_dt <- occ_region_overlap(
          spec.occ = spec_dt,
          region = region,
          region.name = spec.name,
          spec = "species",
          long = "decimalLongitude",
          lat = "decimalLatitude",
          verbose = verbose
        )
      },
      file = nullfile()
    ))

    cols_to_keep <- c("species", "scientificName", "taxonRank", "speciesKey", "totalNobs", "regionNobs", "propNobs")

    unique_counts <- unique(spec_dt, by = "scientificName")
    unique_counts <- unique_counts[, ..cols_to_keep]
    unique_counts[, (paste0("in", region.name)) := fifelse(regionNobs > 0, TRUE, FALSE)]

    spec_dt_out <- rbindlist(list(spec_dt_out, unique_counts), fill = TRUE)

    rm(spec_dt, unique_counts)
    invisible(gc())
  }

  # merge with checklist
  spec_dt_out <- spec_dt_out[checklist, on = c("species", "scientificName")]

  rm(checklist)
  invisible(gc())

  setorder(spec_dt_out, -propNobs)

  catn("Writing file to:", colcat(file.out, color = "output"))

  spec_in_region <- spec_dt_out[regionNobs > 0]
  fwrite(spec_in_region, paste0(dirname(file.out), "/occ-in-region.csv"), bom = TRUE)

  fwrite(spec_dt_out, file.out, bom = TRUE)

  end_timer(time)

  return(spec_dt_out)
}

#' Find the most frequently used name in a vector
#' 
#' @param x A vector containing names to analyze
#' @param max.number Maximum number of unique names allowed before selecting most frequent
#' @param verbose Boolean to control detailed output printing
#' @return The most frequently used name (or unique name if only one exists)
most_used_name <- function(x, col.name, max.number = 1, suppress = FALSE, verbose = FALSE) {
  # Get unique values including empty strings
  name <- unique(x)
  n_unique <- length(name)
  total_records <- length(x)
  
  if (n_unique > max.number) {
    # Print column header for clarity
    vebcat(sprintf("\nFound too many different names for: %s (%d)", col.name, n_unique), veb = !suppress)
    
    # Create frequency table including NAs
    freq_table <- table(x, useNA = "always")
    freq_counts <- as.numeric(freq_table)
    names(freq_counts) <- names(freq_table)
    
    if (verbose) {
      catn("\nFrequency table:")
      freq_df <- data.frame(
        Name = names(freq_counts),
        Count = freq_counts,
        Percentage = sprintf("%.1f%%", 100 * freq_counts / total_records)
      )
      freq_df$Name[is.na(freq_df$Name)] <- "<NA>"
      freq_df <- freq_df[order(-freq_df$Count), ]
      print(freq_df, row.names = TRUE)
    }
    
    # Get the most frequent value
    most_freq_idx <- which.max(freq_counts)
    name <- names(freq_counts)[most_freq_idx]
    freq_count <- freq_counts[most_freq_idx]
    
    # Calculate percentage using total records
    percentage <- sprintf("%.1f", 100 * freq_count / total_records)
    
    vebcat(sprintf("Selected name: '%s' (%d occurrences, %s%%)", 
                 ifelse(is.na(name), "", name), 
                 freq_count,
                 percentage), veb = !suppress)
  }
  
  return(ifelse(is.na(name), "", name))
}

#' Condense taxonomic information from multiple columns
#' 
#' @param spec.dt Data table containing taxonomic columns
#' @param verbose Boolean to control detailed output printing
#' @return Condensed data table with single row of most frequent taxonomic names
condense_taxons <- function(spec.dt, suppress = FALSE, verbose = FALSE) {
  vebcat("Condensing taxon columns.", veb = verbose)
  
  cols_to_select <- c("kingdom", "phylum", "class", "order", "family", 
                      "genus", "species", "infraspecificEpithet", 
                      "taxonRank", "scientificName")
  
  taxon_cols <- spec.dt[, ..cols_to_select]
  ct_cols <- data.table(matrix(ncol = length(cols_to_select), nrow = 1))
  setnames(ct_cols, names(taxon_cols))
  
  cols_to_ignore <- c("infraspecificEpithet", "taxonRank", "scientificName")
  
  for (i in 1:length(taxon_cols)) {
    column <- taxon_cols[[i]]
    col_name <- names(taxon_cols)[i]
    
    # Only show detailed frequency table for non-ignored columns if verbose
    should_show_details <- !(col_name %in% cols_to_ignore)
    
    name <- most_used_name(column,
                           col_name,
                           max.number = 1, 
                           suppress = suppress,
                           verbose = verbose && should_show_details)
    
    ct_cols[[1, i]] <- name
  }
  
  return(ct_cols)
}

handle_wgsrpd_names <- function(wgsrpd.dir = "./resources/region/wgsrpd") {
  vebcat("Handling WGSRPD shapefiles", color = "funInit")

  dirs <- list.dirs(wgsrpd.dir)[-1]

  lvls <- sub(paste0(wgsrpd.dir, "/"), "", dirs)

  for (i in 1:length(dirs)) {
    dir <- dirs[i]
    lvl <- lvls[i]
    shapefile <- paste0(dir, "/", lvl, ".shp")

    catn("Handling WGSRPD", lvl)

    wgsrpd <- load_region(shapefile)

    if (lvl == "level1") {
      names <- c("level1Code", "level1Name")
    } else if (lvl == "level2") {
      names <- c("level2Code", "level1Code", "level1Name", "level2Name")
    } else if (lvl == "level3") {
      names <- c("level3Name", "level3Code", "level2Code", "level1Code")
    } else if (lvl == "level4") {
      names <- c("isoCode", "level4Name", "level4Code", "level4-2", "level3Code", "level2Code", "level1Code")
    }

    names(wgsrpd) <- names

    file.remove(shapefile)

    writeVector(wgsrpd, shapefile, overwrite = TRUE)
  }

  vebcat("WGSRPD shapefiles handled successfully", color = "funSuccess")

  return(invisible())
}

# This function identifies the polygon where the centroid is closest to the species total mean longitude and latitude values. Meaning not the centroid longitude latitude values themselves, but the species total mean coordinates in a given region.
find_wgsrpd_region <- function(spec.dt, projection = "longlat", longitude = "decimalLongitude", latitude = "decimalLatitude", wgsrpd.dir = "./resources/region/wgsrpd", wgsrpdlvl = "3", wgsrpdlvl.name = TRUE, unique = TRUE, suppress = FALSE, verbose = FALSE) {
  vebcat("Identifying WGSRPD regions", color = "funInit", veb = !suppress)

  wgsrpdlvl.shape <- paste0(wgsrpd.dir, "/level", wgsrpdlvl, "/level", wgsrpdlvl, ".shp")

  wgsrpd_lvl <- load_region(wgsrpdlvl.shape)
  
  if (wgsrpdlvl.name) {
    wgsrpdlvl_name <- paste0("level", wgsrpdlvl, "Name")
  } else {
    wgsrpdlvl_name <- paste0("level", wgsrpdlvl, "Code")
  }

  wgsrpd_length <- length(names(wgsrpd_lvl))

  dt <- copy(spec.dt)

  vebprint(names(dt), verbose, "Species data table column names:")

  dt <- dt[!is.na(dt[[longitude]]) & !is.na(dt[[latitude]])]

  vebprint(nrow(dt), verbose, "Species data table row length:")

  prj <- get_crs_config(projection, verbose)
  
  # First get points in region
  vebcat("Converting coordinates to points", veb = !suppress)
  occ_points <- vect(dt, geom = c(longitude, latitude), crs = prj)
  
  # Debugging output
  vebprint(occ_points, verbose, "Occ points before intersect:")
  vebprint(wgsrpd_lvl, verbose, "WGSRPD level data:")

  vebcat("Finding the intersecting points and regions", veb = !suppress)
  wgsrpd_region <- terra::intersect(wgsrpd_lvl, occ_points)
  
  if (nrow(wgsrpd_region) == 0) {
    vebcat("No WGSRPD region intersections found.", veb = !suppress)
    
    # Create a data.table with original coordinates and "unknown" for the region name
    wgsrpd_region_dt <- as.data.table(
      setNames(
        list("unknown", NA, NA, NA, dt[[longitude]], dt[[latitude]]),
        c(
          paste0("level", wgsrpdlvl, "Name"),
          paste0("level", wgsrpdlvl, "Code"),
          paste0("level", as.numeric(wgsrpdlvl) - 1, "Code"),
          paste0("level", as.numeric(wgsrpdlvl) - 2, "Code"),
          "speciesCentroidLong",
          "speciesCentroidLat",
          paste0("level", wgsrpdlvl, "CentroidLong"),
          paste0("level", wgsrpdlvl, "CentroidLat")
        )
      )
    )
    
    if (unique) {
      wgsrpd_region_dt <- unique(wgsrpd_region_dt, by = wgsrpdlvl_name)
    }
  } else { 
      centroids <- get_centroid_subregion(wgsrpd_region, wgsrpdlvl_name, centroid.per.subregion = TRUE, verbose = verbose)
      
        index <- match(values(wgsrpd_region)[[wgsrpdlvl_name]], names(centroids))

        # Here it sets the intersected centroid
        wgsrpd_region[["speciesCentroidLong"]] <- sapply(centroids[index], function(x) terra::crds(x)[1])
        wgsrpd_region[["speciesCentroidLat"]] <- sapply(centroids[index], function(x) terra::crds(x)[2])

    # centroid_dt <- data.table(
    #   level3Name = names(centroids),
    #   speciesCentroidLong = sapply(centroids, \(x) terra::crds(x)[1]),
    #   speciesCentroidLat = sapply(centroids, \(x) terra::crds(x)[2])
    # )
    # 
    # wgsrpd_region_dt <- centroid_dt[wgsrpd_region_dt, on = "level3Name"]
    
    
    region_centroids <- get_centroid_subregion(wgsrpd_lvl, wgsrpdlvl_name, centroid.per.subregion = TRUE, verbose = verbose)
    index_region <- match(values(wgsrpd_region)[[wgsrpdlvl_name]], names(region_centroids))

    # Here it sets the original region centroid
    wgsrpd_region[[paste0("level", wgsrpdlvl, "CentroidLong")]] <- sapply(region_centroids[index_region], function(x) terra::crds(x)[1])
    wgsrpd_region[[paste0("level", wgsrpdlvl, "CentroidLat")]] <- sapply(region_centroids[index_region], function(x) terra::crds(x)[2])
    
    # region_centroid_dt <- data.table(
    #   level3Name = names(region_centroids),
    #   CentroidLong = sapply(region_centroids, function(x) terra::crds(x)[1]),
    #   CentroidLat = sapply(region_centroids, function(x) terra::crds(x)[2])
    # )
    # 
    # # Rename columns dynamically
    # setnames(region_centroid_dt, c("CentroidLong", "CentroidLat"), 
    #          c(paste0("level", wgsrpdlvl, "CentroidLong"), paste0("level", wgsrpdlvl, "CentroidLat")))
    # 
    # wgsrpd_region_dt <- region_centroid_dt[wgsrpd_region_dt, on = "level3Name"]
    
    vebcat("Finishing up", veb = !suppress)
    wgsrpd_region_dt <- as.data.table(wgsrpd_region)
    wgsrpd_region_dt <- wgsrpd_region_dt[, c(
      names(wgsrpd_region_dt)[1:wgsrpd_length],
      "speciesCentroidLong", "speciesCentroidLat",
      paste0("level", wgsrpdlvl, "CentroidLong"),
      paste0("level", wgsrpdlvl, "CentroidLat")
    ), with = FALSE]
    
    if (unique) {
      wgsrpd_region_dt <- unique(wgsrpd_region_dt, by = wgsrpdlvl_name)
    }
  }

  vebcat("All WGSRPD regions identified successfully", color = "funSuccess", veb = !suppress)

  return(wgsrpd_region_dt)
}

# This function finds the centroid coordinate of the wgsrpd level name
get_wgsrpd_polygon_centroid <- function(wgsrpd_dt, wgsrpd.dir = "./resources/region/wgsrpd", wgsrpdlvl = "3") {
  # Load WGSRPD shapefile
  wgsrpdlvl.shape <- paste0(wgsrpd.dir, "/level", wgsrpdlvl, "/level", wgsrpdlvl, ".shp")
  wgsrpd_lvl <- load_region(wgsrpdlvl.shape)
  
  # Get unique regions from input data
  level_name <- paste0("level", wgsrpdlvl, "Name")
  unique_regions <- unique(wgsrpd_dt[[level_name]])
  
  # Initialize results data.table
  result_dt <- data.table(
    level3Name = character(),
    centroidLong = numeric(),
    centroidLat = numeric()
  )
  
  # Calculate geometric centroid for each region
  for (region in unique_regions) {
    
    # Subset polygon for this region
    region_poly <- wgsrpd_lvl[wgsrpd_lvl[[level_name]] == region]
    
    # Get geometric centroid (force inside = TRUE to ensure point is within polygon)
    region_centroid <- terra::centroids(region_poly, inside = TRUE)
    centroid_coords <- terra::crds(region_centroid)
    
    # Add to results
    result_dt <- rbind(result_dt, data.table(
      level3Name = region,
      centroidLong = centroid_coords[1],
      centroidLat = centroid_coords[2]
    ))
  }
  
  return(result_dt)
}

clean_spec_occ <- function(dt, column = "cleanName", projection, resolution, seed, long = "decimalLongitude", lat = "decimalLatitude", verbose = FALSE) {
  vebcat("Getting Long/Lat values.", veb = verbose)

  dt <- dt[!duplicated(dt, by = c(long, lat, column))]

  vebcat("Getting species without any coordinates.", veb = verbose)

  no_coords_sp <- dt[, .N, by = column][N == 0]

  # Check if any coordinates are left for each species
  if (nrow(no_coords_sp) > 0) {
    catn("No coordinates left for the species. Removing it form the data table.")
    return(NULL)
  }

  vebcat("Cleaning species using coordinateCleaner.", veb = verbose)
  # "flag" sets adds true/false to the corresponding tests
  tryCatch(
    {
      dt <- suppressWarnings(
        clean_coordinates(
          dt,
          species = column,
          value = "clean"
        )
      )
    },
    error = function(e) {
      vebcat("Error when cleaning Coordinates:", e$message, color = "nonFatalError")
      stop(e)
    }
  )

  if (nrow(dt) == 0) {
    catn("All species were removed in the Coordinate cleaning process.")
    return(NULL)
  }

  vebcat("Thinning data.", veb = verbose)

  tryCatch(
    {
      thinned_dt <- thin_occ_data(
        dt,
        long = long,
        lat = lat,
        projection = projection,
        res = resolution,
        seed = seed,
        verbose = verbose
      )
    },
    error = function(e) {
      vebcat("Error when thinning data:", e$message, color = "fatalError")
      stop(e)
    }
  )

  vebcat("Calculating median latitude.", veb = verbose)

  tryCatch(
    {
      thinned_dt[, medianLat := calc_abs_by_col(
        data = thinned_dt,
        column = "decimalLatitude",
        method = "median",
        verbose = verbose
      ), by = "cleanName"]
    },
    error = function(e) {
      vebcat("Error when Calculating median latitude:", e$message, color = "fatalError")
      stop(e)
    }
  )

  return(thinned_dt)
}

prepare_species <- function(dt, process.dir, projection = "longlat", verbose = FALSE) {
  if (!is.data.table(dt)) {
    stop("Input must be a data.table")
  }
  vebcat("Preparing species", color = "funInit")

  projection <- get_crs_config(projection)

  vebprint(head(dt, 3), verbose, text = "Data sample:")

  thinned_dt <- clean_spec_occ(
    dt,
    column = "cleanName",
    projection,
    resolution = config$projection$raster_scale_m,
    seed = config$simulation$seed,
    long = "decimalLongitude",
    lat = "decimalLatitude",
    verbose
  )

  if (is.null(thinned_dt)) {
    return(NULL)
  }

  sp_points <- vect(thinned_dt, geom = c("decimalLongitude", "decimalLatitude"), crs = projection)

  if (is.null(sp_points)) {
    catn("Excluding", spec.name, "from further processing.")
    return(list(
      excluded = TRUE
    ))
  }

  return(sp_points)
}

prepare_environment <- function(sp_points, biovars, verbose = FALSE) {
  # Create an empty matrix to store environmental values
  env_values <- matrix(nrow = nrow(sp_points), ncol = terra::nlyr(biovars))
  colnames(env_values) <- names(biovars)
  # set the crs()
  crs(sp_points) <- crs(biovars[[1]])

  if (identical(terra::crs(sp_points), terra::crs(biovars[[1]]))) {
    vebcat("CRS for bio variables and species are identical", color = "proSuccess", veb = verbose)
  } else {
    vebcat("CRS for bio variables and species are not identical:", color = "nonFatalError", veb = verbose)
  }

  # Extract environmental values for each point
  for (i in 1:terra::nlyr(biovars)) {
    vebcat("Extracting biovariable for:", highcat(names(biovars)[i]), veb = verbose)

    vebprint(biovars[[i]], veb = verbose, text = "Biovar sample:")

    env_values[, i] <- terra::extract(biovars[[i]], sp_points, df = TRUE, ID = FALSE)[, 1]
  }

  vebcat("Cleaning extracted data.", veb = verbose)

  if (any(is.na(env_values)) || any(env_values == "")) {
    vebcat("Some extracted values are NA.", color = "nonFatalError", veb = verbose)

    env_values[env_values == ""] <- NA

    clean_env_matrix <- na.omit(env_values)

    if (any(is.na(clean_env_matrix))) {
      vebcat("Cleaning NA-values failed.", color = "nonFatalError", veb = verbose)
    } else {
      vebcat("Cleaning NA-values was successfull.", color = "proSuccess", veb = verbose)
    }
  } else {
    vebcat("No extracted values are NA.", color = "proSuccess", veb = verbose)

    clean_env_matrix <- env_values
  }

  return(clean_env_matrix)
}
