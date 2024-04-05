source_all("./src/hypervolume/data_acquisition/components")

data_acquisition <- function(shapefile, iteration, show.plot, verbose, warn.file, err.file) {
  vebcat("Initiating data acquisition protocol", color = "funInit")
  
  withCallingHandlers(
    {
      biovars_world <- get_wc_data(
        var = climate.var, 
        res = climate.res, 
        show.plot = show.plot, 
        verbose = verbose
      )
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when getting worldClim data", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when getting worldClim data", iteration = iteration)
  )
  
  vebcat("Scaling biovars", veb = verbose)
  
  withCallingHandlers(
    {
      biovars_world <- scale_biovars(
        biovars_world, 
        verbose = verbose
      )
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when scaling biovars_world", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when scaling biovars_world", iteration = iteration)
  )
  
  vebcat("Acquiring biovars_region.", veb = verbose)
  
  withCallingHandlers(
    {
      biovars_region <- get_region_data(
        biovars_world, 
        shapefile = shapefile, 
        projection = "longlat", 
        verbose = verbose
      )
    },
    warning =  function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when acquiring region data", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when acquiring region data", iteration = iteration)
  )
  
  coord_uncertainty <- calc_coord_uncertainty(
    region = biovars_region,
    unit.out = "m",
    dir.out = "./outputs/hypervolume/data_acquisition/logs",
    verbose = verbose
  )
  
  
  vebcat("Data acquisition protocol completed successfully", color = "funSuccess")
  
  return(list(
    world = biovars_world,
    region = biovars_region
  ))
}
