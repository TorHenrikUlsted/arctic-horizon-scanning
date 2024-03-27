reproject_region <- function(region, projection, line_issue = F, show_plot = F, verbose = T) {
  
  region_name <- strsplit(deparse(substitute(region)), "\\$")[[1]][[2]]
  catn("Reprojecting", region_name)
  
  if (line_issue == T) {
    catn("Attempting to fix line issues.")
    
    catn("Getting extents.")
    ext_east <- terra::ext(ext(region)$xmin, 0, ext(region)$ymin, ext(region)$ymax)
    ext_west <- terra::ext(0.00001, ext(region)$xmax, ext(region)$ymin, ext(region)$ymax)
    
    catn("Cropping in half.")
    vect_east <- terra::crop(region, ext_east)
    vect_west <- terra::crop(region, ext_west)
    
    catn("Reprojecting to longlat.")
    proj_east <- terra::project(vect_east, longlat_crs)
    catn("Plotting left side.")
    plot(proj_east)
    proj_west <- terra::project(vect_west, longlat_crs)
    catn("plotting right side.")
    plot(proj_west)
    
    region_longlat <- rbind(proj_west, proj_east)
    
    if (show_plot == T) plot(region)
    
    return(region_longlat)
  }
  
  
  if (projection == "longlat") {
    
    catn("Choosing", highcat("longlat"), "coordinate system.")
    
    if (grepl("+proj=longlat", crs(region, proj = T), fixed = TRUE) == FALSE) {
      catn("Is not longlat, will be reprojected.")
    } 
    
    print(crs(longlat_crs, proj = T))
    prj <- longlat_crs
    
  } else if (projection == "laea") {
    
    catn("Choosing", highcat("polar"), "coordinate system. \n")
    if (grepl("+proj=laea", crs(region, proj = T), fixed = TRUE) == FALSE) {
      vebcat("Is not polar.", color = "nonFatalError")
    } 
    prj <- laea_crs
    
  } else stop("You can only choose projection 'longlat' or 'laea'.")
  
  
  if (!isTRUE(identical(crs(region), prj))) {
    vebcat("Original CRS not identical to current CSR.", veb = verbose)
    vebcat("Reprojecting", highcat(region_name), "to: ", crs(prj, proj = T), veb = verbose) else next
    
    reproj_region <- terra::project(region, prj)
    
    if (!isTRUE(identical(crs(reproj_region), prj))) {
      vebcat("Reprojection completed successfully", color = "funSuccess")
    }  else  {
      vebcat("Reprojection failed.", color = "nonFatalError")
    }
  } else {
    vebcat("Original CRS identical to current CSR.", veb = verbose)
    vebcat(region_name, "reprojected successfully.", color = "funSuccess")
  }
  
  return(reproj_region)
}