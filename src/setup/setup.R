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

  region <- acquire_region(
    shapefiles = c(
      cavm = "./resources/region/cavm2003/cavm.shp",
      glims = "./resources/region/glims_db/glims_polygons.shp"
      # rgi = "./resources/region/rgi/RGI2000-v7.0-C-05_greenland_periphery.shp",
      # rgi_sj = "./resources/region/rgi/sval_jan/RGI2000-v7.0-C-07_svalbard_jan_mayen.shp"
      # wwfEcoRegion = "./resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp"
    )
  )

  region$cavm <- reproject_region(region$cavm, projection = "longlat", line_issue = T)

  region$glims <- reproject_region(region$glims, projection = "longlat")

  identical(crs(region$cavm), crs(region$glims))

  if (any(!is.valid(region$glims))) {
    cat(red("Some geoms are invalid. \n"))
    cat("Attempting to fix \n")
    valid_glims <- makeValid(region$glims)
    if (any(!is.valid(valid_glims))) stop(red("Failed to fix invalid geoms. \n")) else cat(green("Successfully made all geoms valid. \n"))
  } else {
    cat(green("All geoms are valid. \n"))
  }
  
  cat("Cropping glims to cavm extents. \n")
  glims_cavm <- crop(valid_glims, ext(region$cavm))

  cat("Erasing glims from cavm. \n")
  cavm_noice <- erase(region$cavm, glims_cavm)

  create_dir_if("./outputs/setup/region")

  writeVector(cavm_noice, "./outputs/setup/region/cavm_edited.shp")
  
  cat(cc$lightGreen("Region setup completed successfully. \n"))

  return(cavm_noice)
}
