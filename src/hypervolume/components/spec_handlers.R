count_observations <- function(spec.list, dimensions, verbose = FALSE) {
  catn("Removing species with too few observations.")
  
  removed_species <- data.table(
    species = character(0),
    observations = integer(0),
    logObservations = integer(0),
    dimensions = integer(0),
    removed = logical(0),
    meanLat = numeric(0),
    filename = character(0)
  )
  
  for (i in 1:length(spec.list)) {
    cat("\rChecking file", i, "/", length(spec.list))
    spec <- spec.list[[i]]
    spec_name <- basename(gsub(".csv", "", spec))
    spec_dt <- fread(spec, select = "decimalLatitude")
    
    nobs <- nrow(spec_dt)
    
    removed <- FALSE
    meanLat <- 0
    
    if (log(nobs) <= length(dimensions)) {
      removed <- TRUE
    } else {
      meanLat <- mean(spec_dt$decimalLatitude, na.rm = TRUE)
    }
    
    removed_species <- rbind(removed_species, data.table(
      species = gsub("-", " ", spec_name),
      observations = nobs,
      logObservations = round(log(nobs), digits = 3),
      dimensions = length(dimensions),
      removed = removed,
      meanLat = meanLat,
      filename = spec
    ))
  }; catn()
  
  return(removed_species)
}

most_used_name <- function(x, max.number = 1) {
  name <- unique(x)
  
  if (length(name) > max.number) {
    vebcat("Found too many different strings.", color = "nonFatalError")
    
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
    meanLat = mean(decimalLatitude, na.rm = TRUE)), 
    by = "countryCode"]
  
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

prepare_species <- function(dt, process.dir, projection = "longlat", verbose = T) {
  if (!is.data.table(dt)) {
    stop("Input must be a data.table")
  }
  
  prep_dir <- paste0(process.dir, "/species-prep")
  
  create_dir_if(prep_dir)
  
  vebcat("Preparing species", color = "funInit")
  
  vebprint(head(dt, 3), verbose, text = "Data frame sample:")
  
  vebcat("Getting Long/Lat values.", veb = verbose)
  
  dt <- dt[!duplicated(dt, by = c("decimalLongitude", "decimalLatitude", "cleanName"))]
  
  vebcat("Getting species without any coordinates.", veb = verbose)
  
  no_coords_sp <- dt[, .N, by = cleanName][N == 0]
  
  # Check if any coordinates are left for each species
  if(nrow(no_coords_sp) > 0) {
    catn("No coordinates left for the species. Removing it form the data table.")
    return(NULL)
  }
  
  vebcat("Cleaning species using coordinateCleaner.", veb = verbose)
  
  # "flag" sets adds true/false to the corresponding tests
  tryCatch({
    dt <- suppressWarnings( clean_coordinates(dt, species = "cleanName", value = "clean") )
  }, error = function(e) {
    vebcat("Error when cleaning Coordinates:", e$message, color = "nonFatalError")
  })
  
  if (nrow(dt) == 0) {
    catn("All species were removed in the Coordinate cleaning process.")
    return(NULL)
  }
  
  vebcat("Thining coordinates with spThin.", veb = verbose)
  vebcat("max.files written:", length(unique(dt$cleanName)), veb = verbose)
  
  mem_lim <- get_mem_usage("total", "gb") * 0.07
  
  df_size <- object.size(dt)
  
  number_bytes <- 8
  boolean_bytes <- 4
  df_est <- df_size * 2
  small_df_est <- (number_bytes * 2) * 48662
  dist_matrix_est <- number_bytes * nrow(dt) * nrow(dt)
  dist_vect_est <- number_bytes * nrow(dt)
  dist_bool_est <- boolean_bytes * nrow(dt)
  df_out_est <- (number_bytes * 2) * nrow(dt)
  total_est <- (df_est + dist_matrix_est + dist_vect_est + dist_bool_est + df_out_est)
  
  etr <- total_est / 1024^3 * 4
  
  splits <- ceiling(etr / mem_lim)
  
  vebcat("Estimated RAM usage:", etr, veb = verbose)
  vebcat("Max allowed usage", mem_lim, veb = verbose)
  vebcat("Splits needed", splits, veb = verbose)
  
  tryCatch({
    df_list <- split(dt, factor(sort(rank(row.names(dt)) %% splits)))
    
    vebcat("Actual splits made:", length(df_list), veb = verbose) 
    
    sp_thinned_list <- lapply(df_list, function(sp) {
      catn("max.files:", length(unique(sp$cleanName)))
      catn("out.base:", paste0(unique(gsub(" ", "_", sp$cleanName))))
      
      sp_thinned <- suppressWarnings(thin(
        loc.data = sp,
        lat.col = "decimalLatitude",
        long.col = "decimalLongitude",
        spec.col = "cleanName",
        locs.thinned.list.return = TRUE,
        write.files = FALSE,
        max.files = length(unique(sp$cleanName)),
        out.dir = paste0(prep_dir, "/species"),
        out.base = paste0(unique(gsub(" ", "_", sp$cleanName))),
        log.file = paste0(prep_dir, "/thin-log.txt"),
        thin.par = 0.1,
        reps = 1,
      ))
      
      sp_thinned <- sp_thinned[[1]]
      
      sp_thinned$cleanName <- unique(sp$cleanName)
      
      sp_thinned <- sp_thinned %>% select(cleanName, everything())
      
      sp_thinned
    })
    
    vebcat("Combining thinned lists.", veb = verbose)
    df_thinned <- do.call(rbind, sp_thinned_list)
    
  }, error = function(e) {
    vebcat("Error when thinning lists:", e$message, color = "nonFatalError")
  })
  
  if (nrow(df_thinned) == 0) {
    catn("All species were removed in the thinning process.")
    return(NULL)
  }
  
  vebcat("Converting species to points.", veb = verbose)
  
  if (projection == "longlat") {
    prj = longlat_crs
    
  } else if (projection == "laea") {
    prj = laea_crs
    
  } else {
    vebcat("missing projection, using longlat.", veb = verbose)
    prj = longlat_crs
  }
  
  tryCatch({
    sp_points = vect(df_thinned, geom=c("Longitude", "Latitude"), crs = prj)
  }, error = function(e) {
    vebcat("Error when making species into points:", e$message, color = "nonFatalError")
  })
  
  if (verbose == T) {
    if(!any(is.na(sp_points))) {
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
  #set the crs()
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
    
    env_values[, i] <- terra::extract(biovars[[i]], sp_points, df=TRUE, ID = FALSE)[,1]
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
    
    
  }  else {
    vebcat("No extracted values are NA.", color = "proSuccess", veb = verbose)
    
    clean_env_matrix <- env_values
  } 
  
  return(clean_env_matrix)
}