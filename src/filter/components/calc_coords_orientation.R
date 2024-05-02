check_coords_orientation <- function(geom_input) {
  catn("Checking coordinate orientation")
  # Convert the WKT string to an sf object
  if (!inherits(geom_input, "SpatVector")) {
    spat_obj <- vect(geom_input)  
  } else {
    spat_obj <- geom_input
  }
  
  # Extract the coordinates of the polygon
  coords <- geom(spat_obj)
  coords <- coords[,3:4]
  
  # Calculate the signed area using the shoelace formula (positive for counter-clockwise, negative for clockwise)
  signed_area <- 0.5 * sum((coords[,1][1:(nrow(coords) - 1)] * coords[,2][2:nrow(coords)]) -
                           (coords[,1][2:nrow(coords)] * coords[,2][1:(nrow(coords) - 1)]))
  # signed_area <- sum((coords[,1][1:(nrow(coords) - 1)] * coords[,2][2:nrow(coords)]) - 
  #                      (coords[,1][2:nrow(coords)] * coords[,2][1:(nrow(coords) - 1)])) / 2
  # 
  
  # Check the orientation and print a message
  if (signed_area > 0) {
    vebcat("The orientation is counter-clockwise (GBIF friendly).", color = "indicator")
    orientation <- "counter-clockwise"
  } else if (signed_area < 0) {
    vebcat("The orientation is clockwise.", color = "indicator")
    orientation <- "clockwise"
  } else {
    stop("STOP: The WKT does not have a clear orientation. \n")
  }
  
  return(orientation)
}