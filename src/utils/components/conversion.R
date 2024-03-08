to.vector <- function(input, terminate = TRUE, verbose = FALSE) {
  if (verbose) cat("Converting to vector.\n")
  
  if (!is.vector(input)) {
    input <- unlist(input)
    
    if (!is.vector(input)) {
      cat("Failed to convert to vector.\n")
      if (terminate) stop("Check input.")
    } else if (is.vector(input)) {
      if (verbose) cat("successfully converted to vector.\n")
    } else {
      cat("Unlisting the input is not making it turn into a vector. Check input.\n")
    }
  } else {
    if (verbose) cat("Input is already a vector.\n")
  }
  
  return(input)
}

check_crs <- function(object, projection, projection.method, verbose = FALSE) {
  
  if (verbose) {
    cat("Input object:\n")
    print(object)
    
    cat("Input projection:\n")
    print(crs(projection, proj = TRUE))
    
    cat("Input projection.method:\n")
    print(projection.method)
  }
  
  if (!identical(crs(object, proj = TRUE), crs(projection, proj = TRUE))) {
    cat("Reprojecting", cc$lightSteelBlue(as.character(crs(object, proj = TRUE))), "-->", cc$lightSteelBlue(as.character(crs(projection, proj = TRUE))), "\n")
    object <- project(object, projection, method = projection.method)
  }  
  
  return(object)
}
