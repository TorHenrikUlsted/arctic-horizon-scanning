handle_region <- function(region) {

  region <- na.omit(region)
  
  return(region)
}

handle_region_dt <- function(dt) {
  vebcat("Handling region", color = "proInit")
  
  ## Handle how the region is to be used when it is converted to a data.table in the visualsation sequence
  ## Check example for usage
  
 vebcat("Region handled Successfully", color = "proSuccess")
  
  return(dt)
}

setup_region <- function(verbose = FALSE) {
  # This function is used to setup the region
  # Here you need to edit the names of the region as well as what you want to do with it if anything
  
  vebcat("Initiating Region setup.", color = "funInit")
  region_setup_timer <- start_timer("region-setup-timer")
  
  resource_dir <- "./resources/region"
  
  result_shp <- paste0("./outputs/setup/region/shape-name/edited-shape.shp") # edit the shape-name
  
  create_dir_if(dirname(result_shp))
  
  if (!file.exists(result_shp)) {
    
    shapefile <- paste0(resource_dir, "/edited-shape.shp")
    
      download_if(
      out.file = shapefile,
      download.file.ext = "zip",
      download.direct = "shapeDownloadLink", # Direct download link. Can be NULL
      download.page = "shapeDownloadPage" # Fallback url
    )
    
    region <- load_region(
      shapefile, 
      verbose = verbose
    )
    
    # Reproject to desired projection
    region <- reproject_region(
      region, 
      projection = "longlat", 
    )
    
    
    
    ## All your handling here ##
    
    
    
    if (all(terra::is.valid(region))) {
      vebcat("region shape is valid", color = "proSuccess")
    } else {
      vebcat("region shape is invalid", color = "fatalError")
      stop("ERROR: was not able to make a valid shapefile.")
    }
    
    create_dir_if(dirname(result_shp))
    
    catn("Writing vector to file:", colcat(result_shp, color = "output"))
    writeVector(region, result_shp)
  }
  
  end_timer(region_setup_timer)
  
  vebcat("Region setup completed successfully.", color = "funSuccess")
  
  return(result_shp)
}
