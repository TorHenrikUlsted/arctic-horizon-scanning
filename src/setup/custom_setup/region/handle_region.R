handle_region <- function(region) {

  region <- na.omit(region)
  
  return(region)
}

handle_region_dt <- function(dt) {
  # add east to west order
  we <- c(
    "North Alaska - Yukon Territory",
    "Western Alaska",
    "Central Canada",
    "Hudson Bay - Labrador",
    "Ellesmere-North Greenland",
    "Western Greenland",
    "Eastern Greenland",
    "North Iceland",
    "Jan Mayen",
    "North Fennoscandia",
    "Svalbard",
    "Franz Joseph Land",
    "Kanin-Pechora",
    "Polar Ural - Novaya Zemlya",
    "Yamal - Gydan",
    "Anabar - Olenyek",
    "Taimyr - Severnaya Zemlya",
    "Kharaulakh",
    "Yana - Kolyma",
    "West Chukotka",
    "Wrangel Island",
    "South Chukotka",
    "East Chukotka"
  )
  
  we_dt <- data.table(floristicProvince = we, we = 1:length(we))
  
  dt <- dt[we_dt, on = .(floristicProvince)]
  
  return(dt)
}

setup_region <- function(verbose = FALSE) {
  vebcat("Initiating Region setup.", color = "funInit")
  region_setup_timer <- start_timer("region-setup-timer")
  
  resource_dir <- "./resources/region"
  
  greenland_ice <- paste0(resource_dir, "/promice/basalmelt.nc")
  
  glims_zip <- paste0(resource_dir, "/glims-db/glims-db.zip")
  glims_shp <- paste0(resource_dir, "/glims_db/glims_polygons.shp")
  
  result_shp <- paste0(resource_dir, "/cavm-noice/cavm-noice.shp")
  
  if (!file.exists(result_shp)) {
    download_region_if(
      region.file = greenland_ice,
      download.link = "https://dataverse.geus.dk/api/access/datafile/:persistentId?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B",
      download.page = "https://dataverse.geus.dk/file.xhtml?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B&version=2.0"
    )
    
    greenland <- load_region(
      region.file = greenland_ice, 
      verbose = verbose
    )
    
    download_region_if(
      region.file = glims_zip,
      download.link = NULL,
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
    
    # Remove Lake (19), Lagoon (20), Non-Arctic (21)
    #cavm <- cavm[cavm$VEGPHYS != 19, ]
    #cavm <- cavm[cavm$VEGPHYS != 20, ]
    cavm <- cavm[cavm$VEGPHYS != 21, ]
    
    cavm_desc <- fread("./resources/region/cavm2003-desc.csv")
    
    index <- match(cavm$FLOREG, cavm$FLOREG)
    
    cavm$country <- cavm$country[index]
    cavm$floristicProvince <- cavm$floristicProvince[index]
    
    # s <- vect("resources/region/cavm-noice/cavm-noice.shp")
    # 
    # s <- s[s$VEGPHYS != 19, ]
    # s <- s[s$VEGPHYS != 20, ]
    # s <- s[s$VEGPHYS != 21, ]
    
    # s_desc <- fread("./resources/region/cavm2003-desc.csv")
    # 
    # index <- match(s$FLOREG, s$FLOREG)
    # 
    # s$country <- s$country[index]
    # s$floristicProvince <- s$floristicProvince[index]
    
    #writeVector(s, "./outputs/setup/region/cavm-noice.shp")
    
    
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
    
    out_shp <- "./outputs/setup/region/cavm-noice.shp"
    
    catn("Writing vector to file:", colcat(out_shape, color = "output"))
    writeVector(cavm_noice, out_shape)
    
    file.copy(out_shape, result_shp)
  }
  
  end_timer(region_setup_timer)
  
  vebcat("Region setup completed successfully.", color = "funSuccess")
  
  return(result_shp)
}

handle_country_names <- function() {
  glonaf_regions <- fread("./resources/region/glonaf-region-desc.csv")
  
  # Initialize the longitude and latitude columns
  glonaf_regions[, `:=`(longitude = NA_real_, latitude = NA_real_)]
  
  for (i in seq_len(nrow(glonaf_regions))) {
    cat("\rGetting long lat for country", i, "/", nrow(glonaf_regions))
    
    country <- glonaf_regions$countryFixed[i]
    
    # Get the coordinates for the country
    c <- try(map('world', regions = country, plot = FALSE, exact = TRUE), silent = TRUE)
    
    # If the country is found in the maps package
    if (!inherits(c, "try-error")) {
      glonaf_regions[i, `:=`(longitude = mean(c$range[1:2]), latitude = mean(c$range[3:4]))]
    }
  }; catn()
  
  na_regions <- glonaf_regions[is.na(glonaf_regions$longitude), ]
  
  na_regions <- unique(na_regions$countryFixed)
  
  View(na_regions)
  
  return(glonaf_regions)
}
