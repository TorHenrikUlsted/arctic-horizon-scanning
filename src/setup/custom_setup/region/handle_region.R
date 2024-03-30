handle_region <- function(region) {
  region_desc <- fread("./resources/region/cavm2003-desc.csv")
  
  index <- match(region$FLOREG, region_desc$FLOREG)
  
  region$country <- region_desc$country[index]
  region$floristicProvince <- region_desc$floristicProvince[index]
  
  # Remove the ice sheet
  region <- region[region$FLOREG != 0, ]
  region <- na.omit(region)
  
  return(region)
}

setup_region <- function() {
  region_setup_timer <- start_timer("region-setup-timer")
  
  resource_dir <- "./resources/region"
  
  greenlad_ice <- paste0(resource_dir, "/promice/basalmelt.nc")
  
  download_region_if(
    region.file = greenland_ice,
    download.link = "https://dataverse.geus.dk/api/access/datafile/:persistentId?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B",
    download.page = "https://dataverse.geus.dk/file.xhtml?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B&version=2.0"
  )
  
  greenland <- load_region(
    region.file = greenland_ice, 
    verbose = verbose
  )
  
  glims_zip <- paste0(resource_dir, "/glims-db/glims-db.zip")
  glims_shp <- paste0(resource_dir, "/glims_db/glims_polygons.shp")
  
  download_region_if(
    region.file = glims_zip,
    download.link = "https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0272_GLIMS_v1/NSIDC-0272_glims_db_north_20230725_v01.0.zip",
    download.page = "https://nsidc.org/data/nsidc-0272/versions/1"
  )
  
  glims <- load_region(
    glims_shp, 
    verbose = verbose
  )
  
  cavm_shp <- paste0(resource_dir, "/cavm2003/cavm.shp")
  
  cavm <- load_region(
    cavm_shp, 
    verbose = verbose
  )
  
  cavm <- reproject_region(cavm, projection = "longlat", line_issue = TRUE)
  
  catn("Setting up Greenland shape.")
  crs(greenland) <- stere_crs
  
  greenland <- reproject_region(greenland, projection = "longlat")
  
  if (identical(crs(cavm, proj = TRUE), crs(greenland, proj = TRUE))) {
    vebcat("Cavm and Greenland crs are identical.", color = "proSuccess")
  } else {
    vebcat("Cavm and Greenland crs are not identical.", color = "nonFatalError")  } 
  
  green_shape <- as.polygons(greenland[[9]])
  
  valid_gl <- fix_shape(green_shape)
  
  catn("Cropping Greenland Ice to cavm extents.")
  green_cavm <- crop(valid_gl, ext(cavm))
  
  catn("Erasing Greenland Ice from cavm.")
  cavm_noice <- erase(cavm, green_cavm)
  
  catn("Setting up GLIMS database shape.")
  glims <- reproject_region(glims, projection = "longlat")
  
  if (identical(crs(cavm_noice, proj = TRUE), crs(glims, proj = TRUE))) {
    vebcat("Cavm and glims crs are identical.", color = "proSuccess")
  } else {
    vebcat("Cavm and glims crs are not identical.", color = "nonFatalError")
  }
  
  valid_glims <- fix_shape(glims)
  
  catn("Cropping glims to cavm extents.")
  glims_cavm <- crop(valid_glims, ext(cavm_noice))
  
  catn("Erasing glims from cavm.")
  cavm_noice <- erase(cavm_noice, glims_cavm)
  
  create_dir_if("./outputs/setup/region")
  
  out_shape <- "./outputs/setup/region/cavm-noice.shp"
  
  catn("Writing vector to file:", colcat(out_shape, color = "output"))
  writeVector(cavm_noice, out_shape)
  
  end_timer(region_setup_timer)
  
  vebcat("Region setup completed successfully.", color = "funSuccess")
  
  return(cavm_noice)
}
