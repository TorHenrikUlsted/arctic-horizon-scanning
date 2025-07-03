wkt_anti_clockwise <- function(extents) {
  anticlockwise_wkt <- sprintf(
    "POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s))",
    formatC(extents[1], format = "f", digits = 5), formatC(extents[3], format = "f", digits = 5),
    formatC(extents[2], format = "f", digits = 5), formatC(extents[3], format = "f", digits = 5),
    formatC(extents[2], format = "f", digits = 5), formatC(extents[4], format = "f", digits = 5),
    formatC(extents[1], format = "f", digits = 5), formatC(extents[4], format = "f", digits = 5),
    formatC(extents[1], format = "f", digits = 5), formatC(extents[3], format = "f", digits = 5)
  )
  
  # Split the string into individual numbers
  numbers <- strsplit(anticlockwise_wkt, " ")[[1]]
  
  # Remove trailing zeros from each number
  numbers <- sapply(numbers, function(x) {
    x <- gsub("0+$", "", x)
    gsub("\\.$", "", x)
  })
  
  # Combine the numbers back into a single string
  anticlockwise_wkt <- paste(numbers, collapse = " ")
  
  return(anticlockwise_wkt)
}

vect_to_wkt <- function(vect, out.file, min.x = NULL, max.x = NULL, min.y = NULL, max.y = NULL) {
  vebcat("Converting vector to Well-known-text format.", color = "funInit")
  
  if (!inherits(vect, "SpatVector")) {
    stop("The input vect must be a SpatVector.")
  }
  
  xy_params <- list(min.x = min.x, max.x = max.x, min.y = min.y, max.y = max.y)
  
  if (any(!sapply(xy_params, function(x) is.null(x) || is.character(x)))) {
    stop("x and y must be characters or NULL.")
  }
  
  if (is.na(crs(vect)) || crs(vect) == "") {
    catn("No projection found, adding longlat.")
    crs(vect) <- "+proj=longlat +datum=WGS84 +no_defs"
  }
  
  ext_region <- as.vector(ext(vect))
  
  for (i in seq_along(xy_params)) {
    param <- xy_params[[i]]
    param_name <- names(xy_params)[i]
    
    if (!is.null(param)) {
      if (param == "min" && param_name == "min.x") {
        ext_region[1] <- -180
      } else if (param == "max" && param_name == "max.x") {
        ext_region[2] <- 180
      } else if (param == "min" && param_name == "min.y") {
        ext_region[3] <- -90
      } else if (param == "max" && param_name == "max.y") {
        ext_region[4] <- 90
      }
    }
  }
  
  catn("Making GBIF friendly WKT.")
  
  anticlockwise_wkt <- wkt_anti_clockwise(ext_region)
  
  catn("WKT: \n", highcat(anticlockwise_wkt))
  
  out_file <- paste0(out.file, "/", deparse(substitute(vect)), "-wkt.txt")
  
  catn("Writing file to:", colcat(out.file))
  
  write(anticlockwise_wkt, out_file)
  
  vebcat("Vector converted to WKT.", color = "funSuccess")
  
  return(anticlockwise_wkt)
}

combine_WKTs = function(regions, out.file, min_x = F, max_x = F, min_y = F, max_y = F) {
  catn("Combining WKTs")
  
  combined_extents = NULL
  
  # Get extent of each region
  for (region_name in names(regions)) {
    region = regions[[region_name]]
    ext_region = ext(region)
    combined_extents = rbind(combined_extents, ext_region)
  }
  
  # Check if any min or max values are true
  if (min_x) min_x = -180 else min_x = min(combined_extents[,1])
  if (max_x) max_x = 180 else max_x = max(combined_extents[,2])
  if (min_y) min_y = -90 else min_y = min(combined_extents[,3])
  if (max_y) max_y = 90 else max_y = max(combined_extents[,4])
  
  
  # Create GBIF friendly combined polygon
  combined_WKT = wkt_anti_clockwise(combined_extents)
  
  catn("Writing file to:", colcat(out.file, color = "output"))
  
  write(combined_WKT, out.file)
  
  catn("Combined WKT:", combined_WKT)
  
  return(combined_WKT)
}