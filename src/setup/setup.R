source_all("./src/setup/components")

setup_sp <- function(test = T, big_test = F) {
  if (test == T) {
    sp_df <- run_test(big_test)
  } else {
    cat(blue("Initiating full run. \n"))
    sp_df <- run_full()
  }

  return(sp_df)
} # end setup

setup_region <- function() {
  # https://nsidc.org/data/glims/data
  # https://dataverse.geus.dk/dataset.xhtml?persistentId=doi:10.22008/FK2/PLNUEO

  region <- acquire_region(
    shapefiles = c(
      cavm = "./resources/region/cavm2003/cavm.shp",
      cavm_noice = "./resources/region/cavm_e/cavm_edited.shp",
      greenland = "./resources/region/promice/basalmelt_orig.nc"
      #glims = "./resources/region/glims_db/glims_polygons.shp"
    )
  )

  region$cavm <- reproject_region(region$cavm, projection = "longlat", line_issue = T)
  
  crs(region$greenland) <- crs(stere_crs)
  
  region$greenland <- reproject_region(region$greenland, projection = "longlat")
  
  if (identical(crs(region$cavm), crs(region$greenland))) cat(green("Cavm and Greenland crs are identical. \n")) else cat(red("Cavm and Greenland crs are not identical. \n"))
  
  green_shape <- as.polygons(region$greenland[[9]])

  if (any(!is.valid(green_shape))) {
    cat(red("Some geoms of greenland are invalid. \n"))
    cat("Attempting to fix \n")
    valid_gl <- makeValid(green_shape)
    if (any(!is.valid(valid_gl))) stop(red("Failed to fix invalid geoms. \n")) else cat(green("Successfully made all geoms valid. \n"))
  } else {
    cat(green("All greenland geoms are valid. \n"))
    valid_gl <- green_shape
  }
  
  cat("Cropping Greenland Ice to cavm extents. \n")
  green_cavm <- crop(valid_gl, ext(region$cavm))
  
  cat("Erasing Greenland Ice from cavm. \n")
  cavm_noice <- erase(region$cavm_noice, green_cavm)
  
  #region$glims <- reproject_region(region$glims, projection = "longlat")
  
  #if (identical(crs(region$cavm), crs(region$glims))) cat(green("Cavm and glims crs are identical. \n")) else cat(red("Cavm and glims crs are not identical. \n"))

  # if (any(!is.valid(region$glims))) {
  #   cat(red("Some geoms of glims are invalid. \n"))
  #   cat("Attempting to fix \n")
  #   valid_glims <- makeValid(region$glims)
  #   if (any(!is.valid(valid_glims))) stop(red("Failed to fix invalid geoms. \n")) else cat(green("Successfully made all geoms valid. \n"))
  # } else {
  #   cat(green("All glims geoms are valid. \n"))
  # }
  
  # cat("Cropping glims to cavm extents. \n")
  # glims_cavm <- crop(valid_glims, ext(region$cavm))
  # 
  # cat("Erasing glims from cavm. \n")
  # cavm_noice <- erase(region$cavm, glims_cavm)

  create_dir_if("./outputs/setup/region")

  writeVector(cavm_noice, "./outputs/setup/region/cavm_edited.shp")
  
  cat(cc$lightGreen("Region setup completed successfully. \n"))

  return(cavm_noice)
}
