prepare_species <- function(dt, projection, verbose = T) {
  if (!is.data.table(dt)) {
    stop("Input must be a data.table")
  }
  
  vebcat("Preparing species", color = "funInit")
  
  vebprint(head(dt, 3), text = "Data frame sample:")
  
  vebcat("Getting Long/Lat values.", veb = verbose)
  
  dt <- dt[!duplicated(dt, by = .("decimalLongitude", "decimalLatitude", "species"))]
  
  no_coords_sp <- dt[, .N, by = species][N == 0]
  
  # Check if any coordinates are left for each species
  if(nrow(no_coords_sp) > 0) {
    catn("No coordinates left for the species. Removing it form the data table.")
    return(NULL)
  }
  
  vebcat("Cleaning species using coordinateCleaner.", veb = verbose)
  
  prep_dir <- "./outputs/hypervolume/prep"
  
  create_dir_if(prep_dir)
  # "flag" sets adds true/false to the corresponding tests
  tryCatch({
    dt <- suppressWarnings( clean_coordinates(dt, value = "clean") )
  }, error = function(e) {
    vebcat("Error when cleaning Coordinates:", e$message, color = "nonFatalError")
  })
  
  if (nrow(dt) == 0) {
    catn("All species were removed in the Coordinate cleaning process.")
    return(NULL)
  }
  
  vebcat("Thining coordinates with spThin.", veb = verbose)
  vebcat("max.files written:", length(unique(dt$species)), veb = verbose)
  
  mem_lim <- get_mem_usage("total", "gb") * 0.07
  
  str_b <- 56 + (nchar(dt$species[1]) * 4)
  num_b <- 8
  row_b <- str_b + (num_b * 2)
  
  n <- (2.8285*nrow(dt))^2 + nrow(dt) * row_b
  
  etr <- n / 1024^3
  
  splits <- ceiling(etr / mem_lim)
  
  vebcat("Estimated RAM usage:", etr, veb = verbose)
  vebcat("Max allowed usage", mem_lim, veb = verbose)
  vebcat("Splits needed", splits, veb = verbose)
  
  tryCatch({
    df_list <- split(dt, factor(sort(rank(row.names(dt)) %% splits)))
    
    vebcat("Actual splits made:", length(df_list), veb = verbose) 
    
    sp_thinned_list <- lapply(df_list, function(sp) {
      catn("max.files:", length(unique(sp$species)))
      catn("out.base:", paste0(unique(gsub(" ", "_", sp$species))))
      
      sp_thinned <- suppressWarnings(thin(
        loc.data = sp,
        lat.col = "decimalLatitude",
        long.col = "decimalLongitude",
        spec.col = "species",
        locs.thinned.list.return = TRUE,
        write.files = FALSE,
        max.files = length(unique(sp$species)),
        out.dir = paste0(prep_dir, "/species"),
        out.base = paste0(unique(gsub(" ", "_", sp$species))),
        log.file = paste0(prep_dir, "/thin-log.txt"),
        thin.par = 0.1,
        reps = 1,
      ))
      
      sp_thinned <- sp_thinned[[1]]
      
      sp_thinned$species <- unique(sp$species)
      
      sp_thinned <- sp_thinned %>% select(species, everything())
      
      sp_thinned
    })
    
    vebcat("Combining thinned lists.", veb = verbose)
    df_thinned <- do.call(rbind, sp_thinned_list)
    
  }, error = function(e) {
    veb("Error when thinning lists:", e$message, color = "nonFatalError")
  })
  
  if (nrow(df_thinned) == 0) {
    catn("All species were removed in the thinning process.")
    return(NULL)
  }
  
  vebcat("Converting species to points.", veb = verbose)
  
  if (projection == "longlat") {
    prj = longlat_crs
    
  } else if (projection == "laea") {
    prj = laea_crs
    
  } else {
    vebcat("missing projection, using longlat.", veb = verbose)
    prj = longlat_crs
  }
  
  tryCatch({
    sp_points = vect(df_thinned, geom=c("Longitude", "Latitude"), crs = prj)
  }, error = function(e) {
    veb("Error when making species into points:", e$message, color = "nonFatalError")
  })
  
  if (verbose == T) {
    if(!any(is.na(sp_points))) {
      vebcat("No values are NA.", color = "proSuccess") 
    } else {
      vebcat("Some values are NA.", color = "nonFatalError")
    }
  }
  
  return(sp_points)
}

prepare_environment <- function(sp_points, biovars, verbose = T) {
  
  # Create an empty matrix to store environmental values
  env_values <- matrix(nrow = nrow(sp_points), ncol = terra::nlyr(biovars))
  colnames(env_values) <- names(biovars)
  #set the crs()
  crs(sp_points) <- crs(biovars[[1]])
  
  if (identical(terra::crs(sp_points), terra::crs(biovars[[1]]))) {
    vebcat("CRS for bio variables and species are identical", color = "proSuccess", veb = verbose) 
  } else {
    vebcat("CRS for bio variables and species are not identical:", color = "nonFatalError", veb = verbose)
  }
  
  # Extract environmental values for each point
  for (i in 1:terra::nlyr(biovars)) {
    vebcat("Extracting biovariable for:", highcat(names(biovars)), veb = verbose)
    
    env_values[, i] <- terra::extract(biovars[[i]], sp_points, df=TRUE, ID = FALSE)[,1]
  }
  
  
  vebcat("Cleaning extracted data.", veb = verbose)
  
  if (any(is.na(env_values)) || any(env_values == "")) {
    vebcat("Some extracted values are NA.", color = "nonFatalError", veb = verbose)
    
    env_values[env_values == ""] <- NA
    
    clean_env_matrix <- na.omit(env_values)
    
    if (any(is.na(clean_env_matrix))) {
      vebcat("Cleaning NA-values failed.", color = "nonFatalError", veb = verbose) 
    } else {
      vebcat("Cleaning NA-values was successfull.", color = "proSuccess", veb = verbose) 
    }
    
    
  }  else {
    vebcat("No extracted values are NA.", color = "proSuccess", veb = verbose)
    
    clean_env_matrix <- env_values
  } 
  
  return(clean_env_matrix)
}