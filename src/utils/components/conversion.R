to.vector <- function(input, terminate = TRUE, verbose = FALSE) {
  vebcat("Converting to vector.", veb = verbose)
  
  if (!is.vector(input)) {
    input <- unlist(input)
    
    if (!is.vector(input)) {
      vebcat("Failed to convert to vector.", color = "nonFatalError")
      if (terminate) stop("Check input.")
    } else if (is.vector(input)) {
      vebcat("successfully converted to vector.", veb = verbose)
    } else {
      cat("Unlisting the input is not making it turn into a vector. Check input.\n")
    }
  } else {
    vebcat("Input is already a vector.", verbose = TRUE)
  }
  
  return(input)
}

check_crs <- function(object, projection, projection.method, verbose = FALSE) {
  
  vebprint(object, verbose, "Input object:")
  vebprint(crs(projection, proj = TRUE), verbose, "Input projection:")
  vebprint(projection.method, verbose, "Input prjection method:")
  
  if (is.null(crs(object)) || is.null(crs(projection))) {
    stop("Either the object or the projection has an empty CRS.")
  }
  
  if (!identical(crs(object, proj = TRUE), crs(projection, proj = TRUE))) {
    catn("Reprojecting", highcat(as.character(crs(object, proj = TRUE))), "-->", highcat(as.character(crs(projection, proj = TRUE))), "\n")
    object <- terra::project(object, projection, method = projection.method)
  }  
  
  return(object)
}
