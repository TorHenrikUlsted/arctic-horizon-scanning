handle_region <- function(region) {

  region <- na.omit(region)
  
  return(region)
}

handle_region_dt <- function(dt) {
  region_names <- list(
    subRegionCode = "FLOREG",
    subRegionName = "floregName",
    subRegionLong = "floregLong",
    subRegionLat = "floregLat"
  )
  
  for (i in seq_along(region_names)) {
    setnames(dt, region_names[[i]], names(region_names)[i])
  }
  
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
  
  we_dt <- data.table(subRegionName = we, westEast = 1:length(we))
  
  dt <- dt[we_dt, on = .(subRegionName)]
  
  return(dt)
}

setup_region <- function(verbose = FALSE) {
  vebcat("Initiating Region setup.", color = "funInit")
  region_setup_timer <- start_timer("region-setup-timer")
  
  resource_dir <- "./resources/region"
  
  result_shp <- paste0("./outputs/setup/region/cavm-noice/cavm-noice.shp")
  create_dir_if(dirname(result_shp))
  
  if (!file.exists(result_shp)) {
    
    cavm_shp <- paste0(resource_dir, "/cavm2003/cavm.shp") # Add download_if
    
    cavm <- load_region(
      cavm_shp, 
      verbose = verbose
    )
    
    # Reproject first becaus of line issue
    cavm <- reproject_region(
      cavm, 
      projection = "longlat", 
      issue.line = TRUE
    )
    
    catn("Removing non-Arctic regions.")
    # Remove Non-Arctic (21) ice Sheet (0), and glaciers (1)
    cavm <- cavm[cavm$VEGPHYS != 21, ]
    cavm <- cavm[cavm$FLOREG != 0, ]
    cavm <- cavm[cavm$LAND != 1, ]
    
    catn("Converting CAVM to PAF floristic regions.")
    # Give space to the new FLOREG regions
    cavm$FLOREG <- ifelse(cavm$FLOREG >= 11 & cavm$FLOREG <= 23, cavm$FLOREG + 1, cavm$FLOREG)
    
    cavm$FLOREG <- ifelse(cavm$FLOREG >= 14 & cavm$FLOREG <= 24, cavm$FLOREG + 1, cavm$FLOREG)
    
    cavm$ID <- 1:nrow(cavm)
    
    # Split Iceland and Jan Mayen
    jan_mayen <- which(cavm$FLOREG == 10 & cavm$LAND == 6)[4] # row index
    
    cavm$FLOREG[jan_mayen] <- 11
    
    # Split Svalbard and Franz Josef Land
    svalbard_jfl <- subset(cavm, cavm$FLOREG == 13)
    
    svalbard_jfl_geom <- geom(svalbard_jfl) # 41 polygons, use coordinates
    
    index <- svalbard_jfl_geom[svalbard_jfl_geom[, 3] >= 40, ] # Get all polygons above the coordinates
    
    franz_josef <- svalbard_jfl[unique(index[,  1])]
    
    cavm$FLOREG[cavm$ID %in% franz_josef$ID] <- 14
    
    # Combine FLOREG values
    # Get indeces and replace them with a given number
    new_floregs <- list(
      north_alaska = list(old = c(1, 3), new = 1),
      hudson_bay = list(old = c(5, 6), new = 5),
      greenland_west = list(old = c(81, 82, 83, 84), new = 81),
      greenland_east = list(old = c(91, 92, 93), new = 91)
    )
    
    for(nf in names(new_floregs)) {
      indices <- which(cavm$FLOREG %in% new_floregs[[nf]]$old)
      cavm$FLOREG[indices] <- new_floregs[[nf]]$new
    }
    
    # Add region descriptions
    cavm_desc <- fread("./resources/region/cavm2003-desc.csv")
    
    index <- match(cavm$FLOREG, cavm_desc$FLOREG)
    
    cavm$country <- cavm_desc$country[index]
    cavm$floregName <- cavm_desc$floristicProvince[index]
    cavm$floregLong <- cavm_desc$floregLong[index]
    cavm$floregLat <- cavm_desc$floregLat[index]
    
    if (all(terra::is.valid(cavm))) {
      vebcat("Cavm shape is valid", color = "proSuccess")
    } else {
      vebcat("Cavm shape is invalid", color = "fatalError")
      stop("ERROR: was not able to make a valid shapefile.")
    }
    
    create_dir_if(dirname(result_shp))
    
    catn("Writing vector to file:", colcat(result_shp, color = "output"))
    writeVector(cavm, result_shp)
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
