##########################
#       Data.table       #
##########################

merge_and_sum <- function(dt1, dt2, sumCol, by, all = TRUE) {
  merged_dt <- merge(dt1, dt2, by = by, all = all)
  merged_dt[[sumCol]] <- rowSums(merged_dt[, c(paste0(sumCol, ".x"), paste0(sumCol, ".y"))], na.rm = TRUE)

  return(merged_dt)
}

find_min_data_col <- function(file_path, verbose = TRUE) {
  catn("Finding column with the least memory allocation needed.")
  # Read the first 10 rows of the file
  dt <- fread(file_path, nrows = 10)

  # Calculate the total character length of each column
  column_lengths <- sapply(dt, function(x) sum(nchar(as.character(x))))

  # Find the column name with the least character length
  least_data_column <- names(column_lengths)[which.min(column_lengths)]

  catn("Column with the least data:", highcat(least_data_column))

  return(least_data_column)
}

set_df_utf8 <- function(df) {
  for (name in names(df)[sapply(dt, is.character)]) {
    df[[name]] <- enc2utf8(df[[name]])
  }

  return(df)
}

# Used to only standardize, can be used with sapply for certain columns in df
standardize_infraEpithet <- function(spec, verbose = FALSE) {
  res <- spec
  
  for (s in names(standard_infraEpithets)) {
    pattern <- paste0("(?i)\\b", gsub("\\.", "", s), "\\.?\\s+(\\w+)\\b")
    
    if (grepl(pattern, res)) {
      epithet <- tolower(gsub(pattern, "\\1", res, perl = TRUE))
      separator <- ifelse(grepl("\\.$", s), "", ". ")
      replacement <- paste0(standard_infraEpithets[[s]], separator, " ", epithet)
      res <- replacement
      
      break
    }
  }
  
  return(res)
}

remove_designations <- function(spec, verbose = FALSE) {
  # Remove ignored designations
  for (d in ignored_designations) {
    pattern <- paste0("\\b", d, "\\b(?:\\s*\\([^)]+\\))?")
    
    spec <- gsub(pattern, "", spec)
  }
  
  return(trimws(spec))
}

remove_infraEpithet <- function(spec, verbose = FALSE) {
  # Remove ignored designations
  
  spec <- remove_designations(spec = spec, verbose = verbose)
  
  for (d in infraEpithet_designations) {
    # Remove the designation and the name that follows it from the species name
    spec <- gsub(paste0("(\\s*\\(?\\s*(?i)", d, "\\.?\\s+)([^\\)]*)\\)?"), "", spec)
  }

  return(spec)
}

# Get the infraspecificEpithet and standardize it
extract_infraEpithet <- function(spec, verbose = FALSE) {
  res <- ""

  for (d in infraEpithet_designations) {
    pattern <- paste0("(?i)\\b", gsub("\\.", "", d), "\\.?\\s+(\\w+)\\b")

    if (grepl(pattern, spec)) {
      match <- regmatches(spec, regexpr(pattern, spec, perl = TRUE))

      if (length(match) > 0) {
        res <- match[[1]]
        break
      }
    }
  }

  vebprint(res, veb = verbose)

  # Standardize the infraspecific epithet
  if (res != "") {
    res <- standardize_infraEpithet(res, verbose = verbose)
  }

  vebprint(res, veb = verbose)

  return(res)
}

order_by_apg <- function(input, by, verbose = FALSE) {
  # Load apg4
  apg <- fread("./resources/taxon/apg4/apg4.txt")
  
  apg_rank <- apg[taxonRank == by]
  
  apg_vect <- apg_rank$scientificName
  
  non_apg <- setdiff(input, apg_vect)
  
  vebcat("non_apg:", veb = verbose)
  vebprint(non_apg, veb = verbose)
  
  input_apg <- input[!input %in% non_apg]
  
  vebcat("input without non_apgs:", veb = verbose)
  vebprint(input_apg, veb = verbose)
  
  vebcat("Ordering apgs", veb = verbose)
  if (is.vector(input)) {
    ordered_apg <- apg_vect[apg_vect %in% input_apg]
    
    # Add the others
    non_apg_order <- c("Lycopodiales", "Selaginellales", "Equisetales", "Osmundales", "Salviniales", "Cyatheales", "Polypodiales", "Ephedrales", "Pinales")
    
    ordered <- c(non_apg_order, ordered_apg)
    print(ordered)
  } else {
    setorderv(input, cols = apg_rank[[by]])
  }
  
  return(ordered)
}

##########################
#        Spatial         #
##########################

extract_ext_to_dt <- function(raster, value = "value", cells = TRUE) {
  rast_extr <- terra::extract(raster, ext(raster), cells = cells)
  rast_dt <- as.data.table(rast_extr)
  names(rast_dt) <- c("cell", value)

  return(rast_dt)
}

convert_spatial_dt <- function(spatial, verbose = FALSE) {
  dt <- as.data.frame(spatial)
  dt <- as.data.table(dt)

  return(dt)
}

reproject_region <- function(region, projection, issue.line = FALSE, issue.threshold = 0.00001, verbose = FALSE) {
  catn("Reprojecting region")
  
  if (issue.line == T) {
    catn("Attempting to fix line issues.")
    
    catn("Getting extents.")
    ext_east <- terra::ext(ext(region)$xmin, 0, ext(region)$ymin, ext(region)$ymax)
    ext_west <- terra::ext(issue.threshold, ext(region)$xmax, ext(region)$ymin, ext(region)$ymax)
    
    catn("Cropping in half.")
    vect_east <- terra::crop(region, ext_east)
    vect_west <- terra::crop(region, ext_west)
    
    catn("Reprojecting to longlat.")
    proj_east <- terra::project(vect_east, longlat_crs)
    catn("Plotting left side.")
    if (verbose) plot(proj_east)
    proj_west <- terra::project(vect_west, longlat_crs)
    catn("plotting right side.")
    if (verbose) plot(proj_west)
    
    region_longlat <- rbind(proj_west, proj_east)
    
    if (verbose) plot(region_longlat)
    
    return(region_longlat)
  }
  
  
  if (projection == "longlat") {
    catn("Choosing", highcat("longlat"), "coordinate system.")
    prj <- longlat_crs
  } else if (projection == "laea") {
    catn("Choosing", highcat("laea"), "coordinate system.")
    prj <- laea_crs
  } else {
    stop("You can only choose projection 'longlat' or 'laea'.")
  }
  
  
  if (!isTRUE(identical(crs(region, proj = TRUE), crs(prj, proj = TRUE)))) {
    vebcat("Original CRS not identical to current CSR.", veb = verbose)
    catn("Reprojecting region to:\n", highcat(crs(prj, proj = T)))
    
    reproj_region <- terra::project(region, prj)
    
    if (!isTRUE(identical(crs(region, proj = TRUE), crs(prj, proj = TRUE)))) {
      vebcat("Reprojection completed successfully", color = "funSuccess")
    } else {
      vebcat("Reprojection failed.", color = "fatalError")
    }
  } else {
    catn("CRS already correct.")
    reproj_region <- region
  }
  
  return(reproj_region)
}

fix_shape <- function(shape, verbose = FALSE) {
  if (any(!is.valid(shape))) {
    catn("Some geoms of", substitute(deparse(shape)), "are invalid.")
    catn("Attempting to fix.")
    valid_shape <- makeValid(shape)
    if (any(!is.valid(valid_shape))) {
      stop("Failed to fix invalid geoms.")
    } else {
      vebcat("Successfully made all geoms valid.", color = "proSuccess")
    }
    return(valid_shape)
  } else {
    catn("All", substitute(deparse(shape)), "geoms are valid.")

    return(shape)
  }
}

find_peaks <- function(data, prominence = 0.1) {
  # Identify local maxima using diff
  peaks <- which(diff(data) > 0 & diff(c(data, 0)) < 0)

  # Filter based on prominence (difference from neighbors)
  filtered_peaks <- peaks[data[peaks] - data[peaks - 1] >= prominence & data[peaks] - data[peaks + 1] >= prominence]

  return(filtered_peaks)
}

calc_lat_res <- function(lat_res, long_res, latitude = 0, unit.out = "km", verbose = FALSE) {
  
  lat_distance <- lat_res * 111.32
  
  long_deg_size <- 111.32 * cos(latitude * pi / 180)
  long_distance <- long_res * long_deg_size
  
  vebprint(lat_distance, text = "Latitude distance:")
  vebprint(long_distance, text = "Longitude distance:")
  
  # Find the highest resolution
  highest_res <- min(lat_distance, long_distance)
  
  vebprint(highest_res, text = "highest Resolution:")
  
  if (unit.out == "km") {
    return(highest_res)
  } else if (unit.out == "m") {
    return(floor(highest_res * 1000))
  }
}

calc_coord_uncertainty <- function(region, unit.out = "km", dir.out, verbose = FALSE) {
  
  catn("Calculating CoordinateUncertainty.")
  
  out_file <- paste0(dir.out, "/coordinateUncertainty-", unit.out, ".txt")
  
  create_dir_if(dir.out)
  create_file_if(out_file)
  
  if (is.character(region)) {
    region <- rast(region)
  }
  
  res_lat <- terra::res(region)[2]
  res_long <- terra::res(region)[1]
  
  # Get latitude based on northern or southern hemisphere
  region_ext <- terra::ext(region)
  
  vebprint(region_ext, text = "Region Extent:")
  
  if (as.numeric(region_ext[4]) > 0) {
    lat <- as.numeric(region_ext[4]) # northern hemisphere
  } else {
    lat <- as.numeric(region_ext[3]) # Southern hemisphere
  }
  
  vebprint(lat, text = "Latitude:")
  
  max_res <- calc_lat_res(
    res_lat, 
    res_long, 
    lat,
    unit.out = unit.out, 
    verbose = verbose
  )
  
  catn("Writing file to:", colcat(out_file, color = "output"))
  
  writeLines(as.character(max_res), out_file)
  
  catn("Lowest CoordinateUncertainty:", colcat(max_res, color = "indicator"))
  
  return(max_res)
}

load_sp_rast <- function(spec.filename) {
  sp_name <- basename(dirname(spec.filename))

  sp_rast <- terra::rast(spec.filename)
  names(sp_rast) <- sp_name

  return(sp_rast)
}

##########################
#        Objects         #
##########################

get_obj_name <- function(...) {
  sapply(as.list(match.call())[-1], deparse)
}

##########################
#         System         #
##########################

source_all <- function(dir) {
  # Get a list of all .R files in the directory and its subdirectories
  r_files <- list.files(path = dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)

  # Source each file
  lapply(r_files, source)

  cat(length(r_files), "scripts sourced from", dir, "and its subdirectories.\n")
}
