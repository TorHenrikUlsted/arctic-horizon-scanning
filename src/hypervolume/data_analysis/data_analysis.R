source_all("./src/hypervolume/data_analysis/components")

data_analysis <- function(biovars_world, biovars_region, biovars_floreg, hv.dims, method, show.plot, verbose, iteration, warn, err) {
  vebcat("Initiating data analysis protocol", color = "funInit")
  
  if (is.null(hv.dims)) {
    vebcat("hv.dims is NULL, has to be either 'manual', 'auto', or a vector of specific dimensions.", color = "nonFatalError")
  }
  
  withCallingHandlers(
    {
      analyzed_data <- analyze_correlation(biovars_region, file.out = "./outputs/hypervolume/data_analysis/correlation", threshold = 0.5, verbose = T)
    },
    warning = function(w) warn(w, warn_txt = "Warning when analyzing correlation"),
    error = function(e) err(e, err_txt = "Error when analyzing correlation")
  )
  
  withCallingHandlers(
    {
      # choose wanted correlation dimensions
      if (any(hv.dims == "manual")) {
        vbcat("Manual mode", color = "indicator")
        catn("Analyze the correlation matrix in the ./outputs/hypervolume/data_analysis/correlation folder and restart the process with specific dimensions.")
        stop()
      } else if (any(hv.dims == "auto")) {
        stop("Not in use yet, possible to add some functions to calculate the best suited variables/dimensions automatically.")
      } else if (is.vector(hv.dims)) {
        catn("Using dimensions:", as.character(hv.dims))
        hv_dims <- hv.dims
      } else {
        stop("hv.dims is not vector format.\n")
      }
      
      biovars_world <- terra::subset(biovars_world, hv_dims)
      
      biovars_region <- terra::subset(biovars_region, hv_dims)
      
      for (i in seq_along(biovars_floreg)) {
        subset <- biovars_floreg[[i]][[hv_dims]]
        biovars_floreg[[i]] <- subset
      }
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
    biovars_floreg,
    region_hv
  ))
}