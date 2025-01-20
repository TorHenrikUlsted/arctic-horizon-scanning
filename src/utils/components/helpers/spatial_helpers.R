#------------------------#
####      Spatial     ####
#------------------------#

is.spatVector <- function(x) {
  if (inherits(x, "SpatVector")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.spatRaster <- function(x) {
  if (inherits(x, "SpatRaster")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

determine_data_nature <- function(x) {
  dt <- datatype(x)
  
  discrete_datatypes <- c("INT1U", "INT2S", "INT4S", "INT4U", "INT8S", "INT8U", "LOG1S")
  continuous_datatypes <- c("FLT4S", "FLT8S")
  
  if (is.spatVector(x)) {
    return("discrete")  # Always use "near" for SpatVectors
  }
  
  if (all(dt %in% continuous_datatypes)) {
    return("continuous")
  } else if (all(dt %in% discrete_datatypes)) {
    return("discrete")
  } else {
    warning("Mixed or unclear data types. Defaulting to 'near' method.")
    return("discrete")
  }
}

get_crs_config <- function(projection.name, vebose = FALSE) {
  
  try({
    projection <- crs(projection.name)
    return(projection)
  }, silent = TRUE)
  
  
  if (!is.character(projection.name)) {
    catn("Projection.name is not a string. Please chanche the input to a string.")
  }
  # Use config list to dynamically update the input_args if changing config values
  input_args <- config$projection$crs
  
  if (!projection.name %in% names(input_args)) {
    stop(paste("Invalid projection. Valid options are:", paste(names(input_args), collapse = ", ")))
  }
  
  if (inherits(crs(config$projection$crs$longlat, proj=T), "Error")) {
    cat("crs inherit failed")
  }
  
  projection <- input_args[[projection.name]]
  
  catn("Using projection:", projection.name)
  
  return(projection)
}

extract_raster_to_dt <- function(raster, region = NULL, value = "value", cells = TRUE, xy = FALSE, verbose = FALSE) {
  vebcat("Extracting Raster and converting to data table.", veb = verbose)
  
  if (is.null(region)) {
    extract_by <- ext(raster)
  } else {
    extract_by <- region
  }
  
  rast_extr <- terra::extract(raster, extract_by, cells = cells, xy = xy)
  rast_dt <- as.data.frame(rast_extr)
  rast_dt <- as.data.table(rast_extr)
  
  if (is.null(region)) {
    names(rast_dt) <- c("cell", value)
  } else {
    if (xy) {
      names(rast_dt) <- c("ID", value, "cell", "longitude", "latitude")
    } else {
      names(rast_dt) <- c("ID", value, "cell")
    }
  }
  
  return(rast_dt)
}

convert_spatial_dt <- function(spatial, verbose = FALSE) {
  dt <- as.data.frame(spatial)
  dt <- as.data.table(dt)
  
  return(dt)
}

fix_antimeridian <- function(x, threshold = 0.00001, verbose) {
  catn("Fixing antimeridian issue.")
  
  ext_east <- terra::ext(ext(x)$xmin, 0, ext(x)$ymin, ext(x)$ymax)
  ext_west <- terra::ext(threshold, ext(x)$xmax, ext(x)$ymin, ext(x)$ymax)
  
  vect_east <- terra::crop(x, ext_east)
  vect_west <- terra::crop(x, ext_west)
  
  proj_east <- terra::project(vect_east, config$projection$crs$longlat)
  proj_west <- terra::project(vect_west, config$projection$crs$longlat)
  
  if (verbose) {
    plot(proj_east, main = "East side")
    plot(proj_west, main = "West side")
  }
  
  x_fixed <- rbind(proj_west, proj_east)
  
  if (verbose) plot(x_fixed, main = "Combined")
  
  return(x_fixed)
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

thin_occ_data <- function(dt, long = "decimalLongitude", lat = "decimalLatitude", projection = "+proj=longlat +datum=WGS84 +ellps=WGS84", res = 1000, seed = 123, verbose = FALSE) {
  setDT(dt)
  setalloccol(dt)
  
  # Add an ID to the input dt
  dt[, ID := .I]
  
  # get resolution in degrees
  resolution <- res / 111320
  
  # Make into point data
  sp_points <- vect(dt, geom = c(long, lat), crs = projection)
  vebprint(sp_points, verbose, "species points:")
  
  sp_ext <- ext(sp_points)
  
  # Create an empty raster with the same extent as the points
  r_points <- rast(sp_ext, crs = projection, res = resolution)
  values(r_points) <- 1:ncell(r_points)
  vebprint(r_points, verbose, "blank raster with same extent as points:")
  
  if (verbose) {
    plot(r_points)
    plot(sp_points, add = TRUE, col = "red")
  }
  
  # Extract the cell number of a point along with coordinates
  r_dt <- extract_raster_to_dt(r_points, sp_points)
  vebprint(r_dt, verbose, "Raster as data.table:")
  
  # Randomly keep unique points in each cell
  set.seed(seed)
  thinned_data <- r_dt[!is.na(cell), .SD[sample(.N, 1)], by = cell]
  
  # ignore NA cells because it could be points on the border, so do not thin them.
  vebprint(thinned_data, verbose, "Thinned data by choosing a random sample:")
  
  # Merge output with original data table
  thinned_dt <- dt[thinned_data, on = "ID", nomatch = 0]
  
  # Remove unnecessary columns
  thinned_dt[, c("ID", "cell", "value") := NULL]
  
  return(thinned_dt)
}

calc_lat_res <- function(lat_res, long_res, latitude = 0, unit.out = "km", verbose = FALSE) {
  lat_distance <- lat_res * 111.32
  
  long_deg_size <- 111.32 * cos(latitude * pi / 180)
  long_distance <- long_res * long_deg_size
  
  vebprint(lat_distance, text = "Latitude distance:")
  vebprint(long_distance, text = "Longitude distance:")
  
  # Find the highest resolution
  most_precise_res <- min(lat_distance, long_distance)
  
  vebprint(most_precise_res, text = "highest Resolution:")
  
  if (unit.out == "km") {
    return(most_precise_res)
  } else if (unit.out == "m") {
    return(floor(most_precise_res * 1000))
  }
}

calc_coord_uncertainty <- function(region, projection = "longlat", unit.out = "km", dir.out, verbose = FALSE) {
  out_file <- paste0(dir.out, "/coordinateUncertainty-", unit.out, ".txt")
  
  if (file.exists(out_file)) {
    max_res <- as.numeric(readLines(out_file))
  } else {
    catn("Calculating CoordinateUncertainty.")
    
    create_dir_if(dir.out)
    create_file_if(out_file)
    
    if (is.character(region)) {
      region <- rast(region)
    }
    
    if (terra::nlyr(region) > 1) {
      region <- terra::subset(region, 1)
    }
    
    region_ext <- terra::ext(region)
    
    vebprint(region_ext, text = "Region Extent:")
    
    region <- check_crs(region, projection = projection, projection.method = "bilinear")
    
    if (projection == "longlat") {
      res_lat <- terra::res(region)[2]
      res_long <- terra::res(region)[1]
      # Get latitude based on northern or southern hemisphere
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
    } else if (projection == "laea") {
      max_res <- floor(terra::res(region)[1])
      
      if (unit.out == "km") {
        max_res <- (max_res / 1000)
      }
    } else {
      stop("Error: only 'longlat' or 'laea' is available as projection parameters.")
    }
    
    catn("Writing file to:", colcat(out_file, color = "output"))
    
    writeLines(as.character(max_res), out_file)
  }
  
  catn("Lowest CoordinateUncertainty:", colcat(max_res, color = "indicator"))
  
  return(max_res)
}

load_sp_rast <- function(spec.filename) {
  sp_name <- basename(dirname(spec.filename))
  
  sp_rast <- terra::rast(spec.filename)
  names(sp_rast) <- sp_name
  
  return(sp_rast)
}

edit_crs <- function(crs.string, string.key, string.new, verbose = FALSE) {
  string.key <- toupper(string.key)
  
  # Split the keyword if it contains a number
  keyword_parts <- strsplit(string.key, "[[:digit:]]+", perl = TRUE)[[1]]
  keyword_num <- as.numeric(gsub("[^[:digit:]]", "", string.key))
  if (is.na(keyword_num)) keyword_num <- 1
  
  # Split the CRS string into parts
  crs_parts <- strsplit(crs.string, "\n", fixed = TRUE)[[1]]
  
  vebprint(crs_parts, verbose, "CRS parts:")
  
  # Find the parts that start with the keyword
  keyword_parts <- grep(paste0("\\b", keyword_parts, "\\b\\[\""), crs_parts)
  
  vebprint(keyword_parts, verbose, "Keyword parts:")
  
  # Check if the keyword exists in the CRS string
  if (length(keyword_parts) >= keyword_num) {
    # Replace the name following the keyword
    crs_parts[keyword_parts[keyword_num]] <- gsub("(?<=\\[\\\").*?(?=\\\")", string.new, crs_parts[keyword_parts[keyword_num]], perl = TRUE)
    
    # Combine the CRS parts back into a string
    new_crs <- paste(crs_parts, collapse = "\n")
  } else {
    catn("Keyword not found in the string.")
    # If the keyword doesn't exist, return the original CRS string
    new_crs <- crs.string
  }
  
  return(new_crs)
}
# This function gets the polygon closest to the mean of all polygon centroids
get_centroid_subregion <- function(region, region.sub = "subRegion", centroid.per.subregion = FALSE, inside = TRUE, verbose = FALSE) {
  uniq_subregions <- unique(region[[region.sub]])
  
  vebprint(uniq_subregions, verbose, "Unique Sub-Region(s):")
  
  if (centroid.per.subregion) {
    sub_region_centroids <- list()
  } else {
    sub_region_centroids <- vect()
  }
  
  vebcat("Acquiring centroid for ", length(uniq_subregions), "subregion(s)", veb = verbose)
  
  for (i in 1:nrow(uniq_subregions)) {
    sub_region_name <- uniq_subregions[i, ]
    
    vebcat("Acquiring centroid for subregion", sub_region_name, veb = verbose)
    
    vebprint(sub_region_name, verbose, "Sub-Region Name:")
    
    sub_region <- region[region[[region.sub]] == sub_region_name]
    
    vebprint(unique(sub_region[[region.sub]]), verbose, "Actual Subset Region:")
    
    vebprint(sub_region, verbose, "Sub Region:")
    
    all_centroids <- terra::centroids(sub_region, inside = inside)
    
    vebprint(all_centroids, verbose, "All centroids:")
    
    n_centroids <- dim(all_centroids)[1]
    vebprint(n_centroids, verbose, "Dimensions:")
    
    if (n_centroids > 1) {
      all_x <- terra::crds(all_centroids)[, 1]
      all_y <- terra::crds(all_centroids)[, 2]
      
      mean_x <- mean(all_x)
      mean_y <- mean(all_y)
      
      euclidean_distances <- sqrt((all_x - mean_x)^2 + (all_y - mean_y)^2)
      
      centroid <- all_centroids[which.min(euclidean_distances), ]
    } else {
      centroid <- all_centroids
    }
    
    if (!centroid.per.subregion) {
      sub_region_centroids <- rbind(sub_region_centroids, centroid)
    } else {
      sub_region_centroids[[i]] <- centroid
      
      names(sub_region_centroids)[[i]] <- sub_region_name
    }
  }
  
  return(sub_region_centroids)
}