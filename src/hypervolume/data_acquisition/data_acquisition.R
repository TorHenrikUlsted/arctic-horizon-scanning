source_all("./src/hypervolume/data_acquisition/components")

data_acquisition <- function(show.plot, method, verbose, iteration, warn, err) {
  cat(blue("Initiating data acquisition protocol \n"))
  
  if (verbose) cat(blue("Acquiring regions. \n"))
  
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
    if (verbose) cat("Renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(names(cavm_floreg)[i]), "\n")
  }
  
  if (verbose) cat(cc$lightGreen("Region acquisition completed successfully \n"))
  
  if (verbose) cat(blue("Acquiring WorldClim biovariables. \n"))
  
  withCallingHandlers(
    {
      biovars_world <- get_wc_data(show.plot = show.plot, verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when getting worldClim data in iteration"),
    error = function(e) err(e, err_txt = "Error when getting worldClim data in iteration")
  )
  
  if (verbose) cat(cc$lightGreen("Biovars_world acquired successfully \n"))
  
  if (verbose) cat(blue("Scaling biovars_world \n"))
  
  withCallingHandlers(
    {
      biovars_world <- scale_biovars(biovars_world, verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when scaling biovars_world in iteration"),
    error = function(e) err(e, err_txt = "Error when scaling biovars_world in iteration")
  )
  
  if (verbose) cat(cc$lightGreen("Biovars_world scaled successfully \n"))
  
  if (verbose) cat(blue("Acquiring biovars_region \n"))
  
  withCallingHandlers(
    {
      biovars_region <- get_region_data(biovars_world, regions, projection = "longlat", verbose = verbose)
    },
    warning =  function(w) warn(w, warn_txt = "Warning when acquiring region data in iteration"),
    error = function(e) err(e, err_txt = "Error when acquiring region data in iteration")
  )
  
  if (verbose) cat(cc$lightGreen("Biovars_region acquired successfully \n"))
  
  if (verbose) cat(cc$lightGreen("Acquiring biovars for floristic regions \n"))
  
  withCallingHandlers(
    {
      biovars_floreg <- get_region_data(biovars_world, cavm_floreg, projection = "longlat", verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when acquiring floristic region data in iteration"),
    error = function(e) err(e, err_txt = "Error when when acquiring floristic region data in iteration")
  )
  
  if (verbose) cat(cc$lightGreen("Biovars_floreg acquired successfully \n"))
  
  cat(cc$lightGreen("Data acquisition protocol completed successfully \n\n"))
  
  return(list(
    biovars_world,
    biovars_region,
    biovars_floreg
  ))
}