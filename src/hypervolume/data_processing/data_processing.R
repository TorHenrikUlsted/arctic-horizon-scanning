source_all("./src/hypervolume/data_processing/components")

data_processing <- function(sp_df, biovars_world, spec.name, method, points.projection = "longlat", verbose = F, iteration, warn, err) {
  cat(blue("Initiating data processing protocol \n"))
  
  if (verbose) cat(blue("Processing species data.\n"))
  
  withCallingHandlers(
    {
      sp_points <- prepare_species(sp_df, projection = "longlat", verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when preparing species in iteration"),
    error = function(e) err(e, err_txt = "Error when preparing species in iteration")
  )
  
  if (is.null(sp_points)) {
    cat("Excluding", spec.name, "from further processing. \n")
    return(list(excluded = TRUE))
  }
  
  if (verbose) {
    cat("Processed environment data sample: \n")
    print(sp_points)
    cat(cc$lightGreen("Species preperation completed successfully. \n"))
    
    cat(blue("Processing environment data.\n"))
    cat("Using biovars:", cc$lightSteelBlue(names(biovars_world)), "\n")
  }
  
  withCallingHandlers(
    {
      sp_mat <- prepare_environment(sp_points, biovars_world, verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when preparing environment in iteration"),
    error = function(e) err(e, err_txt = "Error when preparing environment in iteration")
  )
  
  if (verbose) {
    cat(cc$lightGreen("Environment preperation completed successfully. \n"))
    cat("Processed environment data sample: \n")
    print(head(sp_mat, 3))
  }
  
  
  cat(cc$lightGreen("Data processing protocol completed successfully. \n\n"))
  
  return(sp_mat)
}