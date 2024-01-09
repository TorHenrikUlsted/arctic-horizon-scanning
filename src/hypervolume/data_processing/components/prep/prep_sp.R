prepare_species <- function(df, projection, verbose = T) {
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or data.table")
  }
  
  if (verbose) cat("Getting Long/Lat values. \n")
  
  df <- df %>% 
    distinct(decimalLongitude, decimalLatitude, species, .keep_all = TRUE)
  
  no_coords_sp <-  df %>% 
    group_by(species) %>% 
    summarise(n = n()) %>% 
    filter(n == 0)
  
  # Check if any coordinates are left for each species
  if(nrow(no_coords_sp) > 0) {
    cat(yellow("Warning: No coordinates left for the species. Removing the species from the data frame. \n"))
    return(NULL)
  }
  
  if (verbose) cat("Cleaning species using coordinateCleaner. \n")
  
  prep_dir <- "./outputs/hypervolume/data-processing/prep"
    
  create_dir_if(prep_dir)
  # "flag" sets adds true/false to the corresponding tests
  df <- suppressWarnings( clean_coordinates(df, value = "clean") )

  if (verbose) cat("Thining coordinates with spThin. \n")

  sp_thinned <- suppressWarnings(thin(
    loc.data = df,
    lat.col = "decimalLatitude",
    long.col = "decimalLongitude",
    spec.col = "species",
    locs.thinned.list.return = TRUE,
    write.files = T,
    max.files = length(unique(df$species)),
    out.dir = paste0(prep_dir, "/species"),
    out.base = paste0(gsub(" ", "_", df$species)),
    log.file = paste(prep_dir, "/thin-log.txt"),
    thin.par = 0.1,
    reps = 1,
  ))
  
  sp_thinned <- sp_thinned[[1]]
  
  sp_thinned$species <- unique(df$species)
  
  sp_thinned <- sp_thinned %>% select(species, everything())
  
  if (verbose == T) cat("Converting species to points \n")
  
  if (projection == "longlat") {
    prj = longlat_crs
    
  } else if (projection == "laea") {
   prj = laea_crs
    
  } else {
    if (verbose) cat(cc$lightCoral("missing projection, using longlat \n"))
    prj = longlat_crs
  }
  
  sp_points = vect(sp_thinned, geom=c("Longitude", "Latitude"), crs = prj)

  if (verbose == T) {
    if(!any(is.na(sp_points))) cat("No values are NA:", green("TRUE") , "\n") else cat("Some values are NA:", red("FALSE") , "\n")   
  }
  
  return(sp_points)
}