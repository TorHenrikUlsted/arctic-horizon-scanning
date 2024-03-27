source_all("./src/hypervolume/data_processing/components")

data_processing <- function(sp_df, biovars_world, spec.name, method, points.projection = "longlat", verbose = F, iteration, warn, err) {
  vebcat("Initiating data processing protocol.", color = "funInit")
  
  vebcat("Starting species data process.", veb = verbose)
  
  withCallingHandlers(
    {
      sp_points <- prepare_species(sp_df, projection = "longlat", verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when preparing species in iteration"),
    error = function(e) err(e, err_txt = "Error when preparing species in iteration")
  )
  
  if (is.null(sp_points)) {
    catn("Excluding", spec.name, "from further processing.")
    return(list(
      excluded = TRUE
    ))
  }
  
  vebprint(sp_points, verbose, "Processed environment data sample:")
  vebcat("Using biovars:", highcat(names(biovars_world)), veb = verbose)
  vebcat("Starting processing of environment data.", veb = verbose)
  
  withCallingHandlers(
    {
      sp_mat <- prepare_environment(sp_points, biovars_world, verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when preparing environment in iteration"),
    error = function(e) err(e, err_txt = "Error when preparing environment in iteration")
  )
  
  vebprint(head(sp_mat, 3), verbose, "Processed environment data sample:")
  
  vebcat("Data processing protocol completed successfully.", color = "funSuccess")
  
  return(sp_mat)
}