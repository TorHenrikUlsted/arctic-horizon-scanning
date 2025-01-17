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


count_observations <- function(spec.list, dimensions, method = "median", verbose = FALSE) {
  catn("Counting species observations.")

  counted_species <- data.table(
    species = character(0),
    observations = integer(0),
    logObservations = integer(0),
    dimensions = integer(0),
    removed = logical(0),
    lat = numeric(0),
    filename = character(0)
  )

  for (i in 1:length(spec.list)) {
    cat("\rChecking file", i, "/", length(spec.list))
    spec <- spec.list[[i]]
    spec_name <- basename(gsub(".csv", "", spec))
    spec_dt <- fread(spec, select = "decimalLatitude")

    nobs <- nrow(spec_dt)

    removed <- FALSE

    if (log(nobs) <= length(dimensions)) {
      removed <- TRUE
    }

    lat <- calc_abs_by_col(spec_dt)

    counted_species <- rbind(counted_species, data.table(
      species = gsub(config$species$file_separator, " ", spec_name),
      observations = nobs,
      logObservations = round(log(nobs), digits = 3),
      dimensions = length(dimensions),
      removed = removed,
      lat = lat,
      filename = spec
    ))
  }
  catn()

  setnames(counted_species, "lat", paste0(method, "Lat"))

  return(counted_species)
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

most_used_name <- function(x, max.number = 1, verbose = FALSE) {
  name <- unique(x)

  if (length(name) > max.number) {
    vebcat("Found too many different names.", color = "nonFatalError")

    freq_table <- table(x)
    vebprint(freq_table, verbose, "Most used names table:")

    name <- names(which.max(freq_table))
  }

  return(name)
}

condense_taxons <- function(spec.dt, verbose = FALSE) {
  catn("Condensing taxon columns.")

  cols_to_select <- c("kingdom", "phylum", "class", "order", "family", "genus", "species", "infraspecificEpithet", "taxonRank", "scientificName")

  taxon_cols <- spec.dt[, ..cols_to_select]
  ct_cols <- data.table(matrix(ncol = length(cols_to_select), nrow = 1)) # condensed taxon cols
  setnames(ct_cols, names(taxon_cols))

  cols_to_ignore <- c("infraspecificEpithet", "taxonRank", "scientificName")

  for (i in 1:length(taxon_cols)) {
    column <- taxon_cols[[i]]

    name <- most_used_name(column, max.number = 1, verbose)

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

# This function identifies the polygon where the centroid is closest to the species total mean longitude and latitude values. Meaning not the centroid longitude latitude values themselves, but the species total mean coordinates.
find_wgsrpd_region <- function(spec.dt, projection = "longlat", longitude = "decimalLongitude", latitude = "decimalLatitude", wgsrpd.dir = "./resources/region/wgsrpd", wgsrpdlvl = "3", wgsrpdlvl.name = TRUE, unique = TRUE, verbose = FALSE) {
  vebcat("Identifying WGSRPD regions", color = "funInit")

  wgsrpdlvl.shape <- paste0(wgsrpd.dir, "/level", wgsrpdlvl, "/level", wgsrpdlvl, ".shp")

  wgsrpd_lvl <- load_region(wgsrpdlvl.shape)

  wgsrpd_length <- length(names(wgsrpd_lvl))

  dt <- copy(spec.dt)

  vebprint(names(dt), verbose, "Species data table column names:")

  dt <- dt[!is.na(dt[[longitude]]) & !is.na(dt[[latitude]])]

  vebprint(nrow(dt), verbose, "Species data table row length:")

  prj <- get_crs_config(projection, verbose)
  
  # First get points in region
  catn("Converting coordinates to points")
  occ_points <- vect(dt, geom = c(longitude, latitude), crs = prj)

  catn("Finding the intersecting points and regions")
  wgsrpd_region <- terra::intersect(wgsrpd_lvl, occ_points)
  
  if (wgsrpdlvl.name) {
    wgsrpdlvl_name <- paste0("level", wgsrpdlvl, "Name")
  } else {
    wgsrpdlvl_name <- paste0("level", wgsrpdlvl, "Code")
  }
  
  # Get the mean location of all points that fall within each WGSRPD region for that species via terra::centroids
  centroids <- get_centroid_subregion(wgsrpd_region, wgsrpdlvl_name, centroid.per.subregion = TRUE, verbose = verbose)
  
  index <- match(values(wgsrpd_region)[[wgsrpdlvl_name]], names(centroids))

  wgsrpd_region[[paste0("level", wgsrpdlvl, "Long")]] <- sapply(centroids[index], function(x) terra::crds(x)[1])

  wgsrpd_region[[paste0("level", wgsrpdlvl, "Lat")]] <- sapply(centroids[index], function(x) terra::crds(x)[2])

  catn("Finishing up")

  wgsrpd_region_dt <- as.data.table(wgsrpd_region)

  wgsrpd_region_dt <- wgsrpd_region_dt[, c(
    names(wgsrpd_region_dt)[1:wgsrpd_length],
    paste0("level", wgsrpdlvl, "Long"),
    paste0("level", wgsrpdlvl, "Lat")
  ), with = FALSE]

  if (unique) {
    wgsrpd_region_dt <- unique(wgsrpd_region_dt, by = wgsrpdlvl_name)
  }

  vebcat("All WGSRPD regions identified successfully", color = "funSuccess")

  return(wgsrpd_region_dt)
}

# This function finds the centroid coordinate of the wgsrpd level name
get_wgsrpd_polygon_centroid <- function(wgsrpd_dt, wgsrpd.dir = "./resources/region/wgsrpd", wgsrpdlvl = "3") {
  # Load WGSRPD shapefile
  wgsrpdlvl.shape <- paste0(wgsrpd.dir, "/level", wgsrpdlvl, "/level", wgsrpdlvl, ".shp")
  wgsrpd_lvl <- load_region(wgsrpdlvl.shape)
  
  # Get unique regions from input data
  unique_regions <- unique(wgsrpd_dt[[paste0("level", wgsrpdlvl, "Name")]])
  
  # Initialize results data.table
  result_dt <- data.table(
    level3Name = character(),
    centroidLong = numeric(),
    centroidLat = numeric()
  )
  
  # Calculate geometric centroid for each region
  for (region in unique_regions) {
    # Subset polygon for this region
    region_poly <- wgsrpd_lvl[wgsrpd_lvl[[paste0("level", wgsrpdlvl, "Name")]] == region]
    
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
