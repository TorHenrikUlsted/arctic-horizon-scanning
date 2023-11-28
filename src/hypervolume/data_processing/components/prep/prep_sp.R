prepare_species <- function(df, projection, verbose = T) {
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or data.table")
  }
  
  if (verbose == T) cat("Getting Long/Lat values. \n")
  
  df <- df %>% 
    distinct(decimalLongitude, decimalLatitude, species, .keep_all = TRUE)
  
  any(duplicated(paste(df$decimalLongitude, df$decimalLatitude)))
  
  no_coords_sp <-  df %>% 
    group_by(species) %>% 
    summarise(n = n()) %>% 
    filter(n == 0)
  
  # Check if any coordinates are left for each species
  if(nrow(no_coords_sp) > 0) {
    cat(yellow("Warning: No coordinates left for the species. Removing the species from the data frame. \n"))
    return(NULL)
  }
  
  cat("Cleaning species using coordinateCleaner. \n")
  
  create_dir_if("./outputs/data_processing/prep")
  # "flag" sets adds true/false to the corresponding tests
  df <- ?clean_coordinates(df, value = "clean")
  
  cat("Thining coordinates with spThin. \n")
  
  sp_thinned <- suppressWarnings(thin(
    loc.data = df,
    lat.col = "decimalLatitude",
    long.col = "decimalLongitude",
    spec.col = "species",
    locs.thinned.list.return = TRUE,
    write.files = T,
    max.files = length(unique(df$species)),
    out.dir = "./outputs/data_processing/prep",
    out.base = paste0(gsub(" ", "_", df$species)),
    log.file = "./outputs/data_processing/prep/spatial_thin_log.txt",
    thin.par = 0.1,
    reps = 1,
  ))
  
  sp_thinned <- sp_thinned[[1]]
  
  sp_thinned$species <- unique(df$species)
  
  sp_thinned <- sp_thinned %>% select(species, everything())
  
  if (verbose == T) cat("Converting species to points \n")
  
  if (projection == "longlat") {
    prj = crs("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84")
    
  } else if (projection == "laea") {
   prj = crs("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
    
  } else {
    cat(red("missing projection, using longlat"))
    prj = crs("+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84")
  }
  
  sp_points = vect(sp_thinned, geom=c("Longitude", "Latitude"), crs = prj)

  if (verbose == T) {
    if(!any(is.na(sp_points))) cat("No values are NA:", green("TRUE") , "\n") else cat("Some values are NA:", red("FALSE") , "\n")   
  }
  
  return(sp_points)
}