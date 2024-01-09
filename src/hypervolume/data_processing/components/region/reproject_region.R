reproject_region <- function(region, projection, line_issue = F, show_plot = F, verbose = T) {
  
  region_name <- strsplit(deparse(substitute(region)), "\\$")[[1]][[2]]
  cat(blue("Reprojecting", region_name, "\n"))
  
  if (line_issue == T) {
    cat("Attempting to fix line issues. \n")
    
    cat("Getting extents. \n")
    ext_east <- terra::ext(ext(region)$xmin, 0, ext(region)$ymin, ext(region)$ymax)
    ext_west <- terra::ext(0.00001, ext(region)$xmax, ext(region)$ymin, ext(region)$ymax)
    
    cat("Cropping in half. \n")
    vect_east <- terra::crop(region, ext_east)
    vect_west <- terra::crop(region, ext_west)
    
    cat("Reprojecting to longlat. \n")
    proj_east <- terra::project(vect_east, longlat_crs)
    cat("Plotting left side \n")
    plot(proj_east)
    proj_west <- terra::project(vect_west, longlat_crs)
    cat("plotting right side \n")
    plot(proj_west)
    
    region_longlat <- rbind(proj_west, proj_east)
    
    if (show_plot == T) plot(region)
    
    return(region_longlat)
  }
  
  
  if (projection == "longlat") {
    
    cat("Choosing", cc$lightSteelBlue("longlat"),  "coordinate system. \n")
    if (grepl("+proj=longlat", crs(region, proj = T), fixed = TRUE) == FALSE) cat(red("Is not longlat, needs reprojection. \n"))
    print(crs(longlat_crs, proj = T))
    prj <- longlat_crs
    
  } else if (projection == "laea") {
    
    cat("Choosing", cc$lightSteelBlue("polar"), "coordinate system. \n")
    if (grepl("+proj=laea", crs(region, proj = T), fixed = TRUE) == FALSE) cat(red("Is not polar \n"))
    prj <- laea_crs
    
  } else stop("You can only choose projection 'longlat' or 'laea'.")
  
  
  if (!isTRUE(identical(crs(region), prj))) {
    if (verbose == T) cat(yellow("Original CRS not identical to current CSR. \n"))
    if (verbose == T) cat("Reprojecting", cc$lightSteelBlue(region_name), "to: ", crs(prj, proj = T), "\n") else next
    
    reproj_region <- terra::project(region, prj)
    
    if (!isTRUE(identical(crs(reproj_region), prj))) cat(green("Reprojection completed successfully: "), crs(reproj_region, proj = T), "\n") else cat(red("Reprojection failed."))
  } else {
    if (verbose == T) cat("Original CRS identical to current CSR. \n")
  }
  
  cat(cc$lightGreen(region_name, "reprojected successfully. \n"))
  
  return(reproj_region)
}