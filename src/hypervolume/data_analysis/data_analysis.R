source_all("./src/hypervolume/data_analysis/components")

data_analysis <- function(biovars_world, biovars_region, hv.dims, method, show.plot, verbose, iteration, warn, err) {
  vebcat("Initiating data analysis protocol", color = "funInit")
  
  if (is.null(hv.dims)) {
    vebcat("hv.dims is NULL, has to be a vector of specific biovariables", color = "nonFatalError")
  }
  
  withCallingHandlers(
    {
      
      biovars_world <- terra::subset(biovars_world, hv_dims)
      
      biovars_region <- terra::subset(biovars_region, hv_dims)
      
    },
    warning = function(w) warn(w, warn_txt = "Warning when subsetting regions"),
    error = function(e) err(e, err_txt = "Error when when subsetting regions")
  )
  
  withCallingHandlers(
    {
      region_hv <- setup_region_hv(biovars_region, out.dir = "./outputs/hypervolume/data_analysis/region", name = "cavm", method = method)
    },
    warning = function(w) warn(w, warn_txt = "Warning when setting up region hypervolume"),
    error = function(e) err(e, err_txt = "Error when setting up region hypervolume")
  )
  
  vebcat("Data analysis protocol completed successfully", color = "funSuccess")
  
  return(list(
    biovars_world,
    biovars_region,
    region_hv
  ))
}