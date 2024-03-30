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

vect_to_wkt <- function(vect, out.file, min.x = F, max.x = F, min.y = F, max.y = F) {
  vebcat("Converting vector to Well-known-text format.", color = "funInit")
  
  if (class(vect) != "SpatVector") {
    stop("The input vect must be a SpatVector.")
  }
  
  if (is.na(crs(vect)) || crs(vect) == "") {
    catn("No projection found, adding longlat.")
    # Assign a CRS if it doesn't have one
    crs(vect) <- "+proj=longlat +datum=WGS84 +no_defs"
  }
  
  ext_region <- ext(vect)
  if (class(ext_region) == "SpatExtent") {
    ext_region <- c(ext_region$xmin, ext_region$xmax, ext_region$ymin, ext_region$ymax)
  }
  
  # Check if any min or max values are true
  if (min.x) ext_region[1] <- -180 
  if (max.x) ext_region[2] <- 180 
  if (min.y) ext_region[3] <- -90 
  if (max.y) ext_region[4] <- 90 
  
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