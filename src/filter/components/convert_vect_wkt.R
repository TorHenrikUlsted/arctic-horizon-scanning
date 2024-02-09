vect_to_wkt <- function(vect, min.x = F, max.x = F, min.y = F, max.y = F) {
  cat(blue("Converting vector to Well-known-text format. \n"))
  
  if (class(vect) != "SpatVector") {
    stop("The input vect must be a SpatVector.")
  }
  
  if (is.na(crs(vect)) || crs(vect) == "") {
    cat("No projection found, adding longlat \n")
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
  
  cat("Making GBIF friendly WKT \n")
  
  anticlockwise_wkt <- sprintf(
    "POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s))",
    formatC(ext_region[1], format = "f", digits = 5), formatC(ext_region[3], format = "f", digits = 5),
    formatC(ext_region[2], format = "f", digits = 5), formatC(ext_region[3], format = "f", digits = 5),
    formatC(ext_region[2], format = "f", digits = 5), formatC(ext_region[4], format = "f", digits = 5),
    formatC(ext_region[1], format = "f", digits = 5), formatC(ext_region[4], format = "f", digits = 5),
    formatC(ext_region[1], format = "f", digits = 5), formatC(ext_region[3], format = "f", digits = 5)
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
  
  cat("WKT: \n", cc$lightSteelBlue(anticlockwise_wkt), "\n")
  
  create_dir_if("./outputs/filter/region")
   
  write(anticlockwise_wkt, paste0("./outputs/filter/region/", deparse(substitute(vect)),"-wkt.txt"))
  
  cat(cc$lightGreen("Region converted to WKT. \n"))
  
  return(anticlockwise_wkt)
}