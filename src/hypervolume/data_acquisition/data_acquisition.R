source_all("./src/hypervolume/data_acquisition/components")

data_acquisition <- function(shapefiles, show.plot, method, verbose, iteration, warn, err) {
  vebcat("Initiating data acquisition protocol", color = "funInit")
  
  vebcat("Loading regions.", veb = verbose)
  
  shapefiles <- c(
    cavm = "./resources/region/cavm-noice/cavm-noice.shp"
  )
  withCallingHandlers(
    {
      regions <- import_regions(shapefiles, "./outputs/hypervolume/data_acquisition/region")
    }, 
    warning = function(w) warn(w, warn_txt = "Warning when importing regions in iteration"),
    error = function(e) err(e, err_txt = "Error when importing regions in iteration")
  )
  
  cavm_floreg <- terra::split(regions$cavm, regions$cavm$FLOREG)
  
  for (i in seq_along(cavm_floreg)) {
    names(cavm_floreg)[i] <- paste("floreg", unique(cavm_floreg[[i]]$FLOREG), sep = "_")
    vebcat("Renaming item", highcat(i), "to", highcat(names(cavm_floreg)[i]), veb = verbose)
  }
  
  vebcat("Loading World Clim data.", veb = verbose)
  
  withCallingHandlers(
    {
      biovars_world <- get_wc_data(show.plot = show.plot, verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when getting worldClim data in iteration"),
    error = function(e) err(e, err_txt = "Error when getting worldClim data in iteration")
  )
  
  vebcat("Scaling biovars", veb = verbose)
  
  withCallingHandlers(
    {
      biovars_world <- scale_biovars(biovars_world, verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when scaling biovars_world in iteration"),
    error = function(e) err(e, err_txt = "Error when scaling biovars_world in iteration")
  )
  
  vebcat("Acquiring biovars_region.", veb = verbose)
  
  withCallingHandlers(
    {
      biovars_region <- get_region_data(biovars_world, regions, projection = "longlat", verbose = verbose)
    },
    warning =  function(w) warn(w, warn_txt = "Warning when acquiring region data in iteration"),
    error = function(e) err(e, err_txt = "Error when acquiring region data in iteration")
  )
  
  vebcat("Acquiring biovars for floristic regions.", veb = verbose)
  
  withCallingHandlers(
    {
      biovars_floreg <- get_region_data(biovars_world, cavm_floreg, projection = "longlat", verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when acquiring floristic region data in iteration"),
    error = function(e) err(e, err_txt = "Error when when acquiring floristic region data in iteration")
  )
  
  vebcat("Data acquisition protocol completed successfully", color = "funSuccess")
  
  return(list(
    biovars_world,
    biovars_region,
    biovars_floreg
  ))
}