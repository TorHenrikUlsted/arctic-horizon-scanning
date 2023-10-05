combine_wkt_anticlockwise <- function(regions, min_x = F, max_x = F, min_y = F, max_y = F) {
  combined_extents <- NULL
  
  cat("Combining the WKTs of: ", names(regions), "\n")
  
  # Get extent of each region
  for (region_name in names(regions)) {
    region <- regions[[region_name]]
    ext_region <- ext(region)
    if (class(ext_region) == "SpatExtent") {
      ext_region <- c(ext_region$xmin, ext_region$xmax, ext_region$ymin, ext_region$ymax)
    }
    combined_extents <- rbind(combined_extents, ext_region)
  }
  
  
  # Check if any min or max values are true
  if (min_x) min_x <- -180 else min_x <- min(combined_extents[, 1])
  if (max_x) max_x <- 180 else max_x <- max(combined_extents[, 2])
  if (min_y) min_y <- -90 else min_y <- min(combined_extents[, 3])
  if (max_y) max_y <- 90 else max_y <- max(combined_extents[, 4])
  
  
  # Create anticlockwise (GBIF friendly) combined polygon
  combined_WKT_anticlockwise <- sprintf(
    "POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s))",
    formatC(min_x, format = "f", digits = 5), formatC(min_y, format = "f", digits = 5),
    formatC(max_x, format = "f", digits = 5), formatC(min_y, format = "f", digits = 5),
    formatC(max_x, format = "f", digits = 5), formatC(max_y, format = "f", digits = 5),
    formatC(min_x, format = "f", digits = 5), formatC(max_y, format = "f", digits = 5),
    formatC(min_x, format = "f", digits = 5), formatC(min_y, format = "f", digits = 5)
  )
  
  # Split the string into individual numbers
  numbers <- strsplit(combined_WKT_anticlockwise, " ")[[1]]
  # Remove trailing zeros from each number
  numbers <- sapply(numbers, function(x) {
    x <- gsub("0+$", "", x)
    gsub("\\.$", "", x)
  })
  # Combine the numbers back into a single string
  combined_WKT_anticlockwise <- paste(numbers, collapse = " ")
  cat("Combined WKT: ", combined_WKT_anticlockwise, "\n")
  write(combined_WKT_anticlockwise, "./outputs/filtering/gbif_retrieval_process/combined_WKT_anticlockwise.txt")
  
  return(combined_WKT_anticlockwise)
}