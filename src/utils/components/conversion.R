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
  tryCatch({
   prj_crs <- crs(projection, proj = TRUE)
  },
    error = function(e) {
      vebcat("Error when trying to get crs.", color = "fatalError")
      catn("Example of Accepted format:")
      catn("+proj=laea +lon_0=0 +lat_0=90 +datum=WGS84")
      catn("Example of not accepted format:")
      catn("laea")
      stop("Check if input projection is a real crs and not a string name.")
    }
  )
  
  vebprint(object, verbose, "Input object:")
  vebprint(crs(projection, proj = TRUE), verbose, "Input projection:")
  vebprint(projection.method, verbose, "Input projection method:")
  
  if (is.null(crs(object)) || is.null(crs(projection))) {
    stop("Either the object or the projection has an empty CRS.")
  }
  
  if (!identical(crs(object, proj = TRUE), prj_crs)) {
    catn("Reprojecting", highcat(as.character(crs(object, proj = TRUE))), "-->", highcat(as.character(crs(projection, proj = TRUE))), "\n")
    object <- terra::project(object, projection, method = projection.method)
  }  
  
  return(object)
}

to_char <- function(input, string = "Updated input:", terminate = TRUE, verbose = FALSE) {
  
  call_stack <- sys.calls()
  is_nested <- length(call_stack) > 1 && !identical(call_stack[[length(call_stack) - 1]][[1]], as.name("do.call"))
  
  tryCatch({
    vebprint(input, verbose, "Original input:")
    
    if (is.character(input)) {
      input.updated <- input
    } else if (is_nested) {
      parent_call <- call_stack[[length(call_stack) - 1]]
      for (i in seq_along(parent_call)) {
        if (identical(eval(parent_call[[i]], parent.frame(2)), input)) {
          input.updated <- deparse(parent_call[[i]])
        }
      }
    } else {
      input.updated <- deparse(substitute(input))
    }
    
    vebprint(input.updated, verbose, string)
    
  }, error = function(e) {
    vebcat("Error: input is not a string. Please input as a string.", color = "fatalError")
    try(vebprint(class(input), text = "Found class:"))
    if (terminate) stop("Edit the input paramter.")
  })
  
  return(input.updated)
}
