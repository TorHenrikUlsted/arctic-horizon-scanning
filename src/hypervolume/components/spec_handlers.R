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

most_used_name <- function(x, max.number = 1) {
  name <- unique(x)

  if (length(name) > max.number) {
    vebcat("Found too many different names.", color = "nonFatalError")

    freq_table <- table(x)

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

    name <- most_used_name(column)

    ct_cols[[1, i]] <- name
  }

  return(ct_cols)
}

condense_country <- function(spec.dt, verbose = FALSE) {
  catn("Condensing country columns.")

  cols_to_select <- c("countryCode", "decimalLongitude", "decimalLatitude")

  country_cols <- spec.dt[, ..cols_to_select]

  country_coords <- country_cols[, .(
    meanLong = mean(decimalLongitude, na.rm = TRUE),
    meanLat = mean(decimalLatitude, na.rm = TRUE)
  ),
  by = "countryCode"
  ]

  code_list <- list(
    XK = "Kosovo"
  )

  country_coords[, country := sapply(countryCode, function(code) {
    if (code %in% names(code_list)) {
      return(code_list[[code]])
    } else if (code %in% countrycode::codelist$iso2c) {
      return(countrycode(code, "iso2c", "country.name"))
    } else if (is.na(code) || is.null(code) || code == "") {
      return(NA)
    } else {
      warning(paste("Invalid countryCode:", code))
      return(NA)
    }
  })]

  # sort
  setcolorder(country_coords, c("countryCode", "country", "meanLong", "meanLat"))

  return(country_coords)
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
  
  catn("Converting coordinates to points")
  occ_points <- vect(dt, geom = c(longitude, latitude), crs = prj)
  
  catn("Finding the intersecting points and regions")
  wgsrpd_region <- terra::intersect(wgsrpd_lvl, occ_points)
  
  if (wgsrpdlvl.name) {
    wgsrpdlvl_name <- paste0("level", wgsrpdlvl, "Name")
  } else {
    wgsrpdlvl_name <- paste0("level", wgsrpdlvl, "Code")
  }
  
  centroids <- get_centroid_subregion(wgsrpd_region, wgsrpdlvl_name, centroid.per.subregion = TRUE, verbose = verbose)
  
  index <- match(values(wgsrpd_region)[[wgsrpdlvl_name]], names(centroids))
  
  wgsrpd_region[[paste0("level", wgsrpdlvl,"Long")]] <- sapply(centroids[index], function(x) terra::crds(x)[1])
  
  wgsrpd_region[[paste0("level", wgsrpdlvl,"Lat")]] <- sapply(centroids[index], function(x) terra::crds(x)[2])
  
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

prepare_species <- function(dt, process.dir, projection = "longlat", verbose = T) {
  if (!is.data.table(dt)) {
    stop("Input must be a data.table")
  }

  vebcat("Preparing species", color = "funInit")

  vebprint(head(dt, 3), verbose, text = "Data frame sample:")

  vebcat("Getting Long/Lat values.", veb = verbose)

  dt <- dt[!duplicated(dt, by = c("decimalLongitude", "decimalLatitude", "cleanName"))]

  vebcat("Getting species without any coordinates.", veb = verbose)

  no_coords_sp <- dt[, .N, by = cleanName][N == 0]

  # Check if any coordinates are left for each species
  if (nrow(no_coords_sp) > 0) {
    catn("No coordinates left for the species. Removing it form the data table.")
    return(NULL)
  }

  vebcat("Cleaning species using coordinateCleaner.", veb = verbose)

  # "flag" sets adds true/false to the corresponding tests
  tryCatch(
    {
      dt <- suppressWarnings(clean_coordinates(dt, species = "cleanName", value = "clean"))
    },
    error = function(e) {
      vebcat("Error when cleaning Coordinates:", e$message, color = "nonFatalError")
    }
  )

  if (nrow(dt) == 0) {
    catn("All species were removed in the Coordinate cleaning process.")
    return(NULL)
  }

  vebcat("Converting species to points.", veb = verbose)

  if (projection == "longlat") {
    prj <- config$projection$crs$longlat
  } else if (projection == "laea") {
    prj <- config$projection$crs$laea
  } else {
    vebcat("missing projection, using longlat.", veb = verbose)
    prj <- config$projection$crs$longlat
  }

  tryCatch(
    {
      thinned_dt <- thin_occ_data(
        dt,
        long = "decimalLongitude",
        lat = "decimalLatitude",
        projection = prj,
        res = config$projection$raster_scale_m,
        seed = config$simulation$seed,
        verbose = verbose
      )

      sp_points <- vect(thinned_dt, geom = c("decimalLongitude", "decimalLatitude"), crs = prj)
    },
    error = function(e) {
      vebcat("Error when thinning data and making species into points:", e$message, color = "nonFatalError")
    }
  )

  if (verbose == T) {
    if (!any(is.na(sp_points))) {
      vebcat("No values are NA.", color = "proSuccess")
    } else {
      vebcat("Some values are NA.", color = "nonFatalError")
    }
  }

  return(sp_points)
}

prepare_environment <- function(sp_points, biovars, verbose = T) {
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
