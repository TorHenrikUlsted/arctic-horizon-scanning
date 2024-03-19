prepare_species <- function(df, projection, verbose = T) {
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or data.table")
  }
  cat("Data frame sample:\n")
  print(head(df, 3))
  
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
  
  prep_dir <- "./outputs/hypervolume/data_processing/prep"
    
  create_dir_if(prep_dir)
  # "flag" sets adds true/false to the corresponding tests
  tryCatch({
    df <- suppressWarnings( clean_coordinates(df, value = "clean") )
  }, error = function(e) {
    cat("Error when cleaning Coordinates:", e$message, "\n")
  })
  
  if (nrow(df) == 0) {
    cat("All species were removed in the Coordinate cleaning process. \n")
    return(NULL)
  }

  if (verbose) cat("Thining coordinates with spThin. \n")
  if (verbose) cat("max.files written:", length(unique(df$species)), "\n")

  mem_lim <- get_mem_usage("total", "gb") * 0.07
  
  str_b <- 56 + (nchar(df$species[1]) * 4)
  num_b <- 8
  row_b <- str_b + (num_b * 2)
  
  n <- (2.8285*nrow(df))^2 + nrow(df) * row_b
  
  etr <- n / 1024^3
  
  splits <- ceiling(etr / mem_lim)
  
  if (verbose) cat("Estimated RAM usage:", etr, "\n")
  if (verbose) cat("Max allowed usage", mem_lim, "\n")
  if (verbose) cat("Splits needed", splits, "\n")
  
  tryCatch({
    df_list <- split(df, factor(sort(rank(row.names(df)) %% splits)))
    
    if (verbose) cat("Actual splits made:", length(df_list), "\n") 
    
    sp_thinned_list <- lapply(df_list, function(sp) {
      cat("max.files:", length(unique(sp$species)), "\n")
      cat("out.base:", paste0(unique(gsub(" ", "_", sp$species))), "\n")

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
    
    if (verbose) cat("Combining thinned lists. \n")
    df_thinned <- do.call(rbind, sp_thinned_list)
    
  }, error = function(e) {
    cat("Error when thinning lists:", e$message, "\n")
  })
  
  if (nrow(df_thinned) == 0) {
    cat("All species were removed in the thinning process. \n")
    return(NULL)
  }
  
  if (verbose == T) cat("Converting species to points \n")
  
  if (projection == "longlat") {
    prj = longlat_crs
    
  } else if (projection == "laea") {
   prj = laea_crs
    
  } else {
    if (verbose) cat(cc$lightCoral("missing projection, using longlat \n"))
    prj = longlat_crs
  }
  
  tryCatch({
    sp_points = vect(df_thinned, geom=c("Longitude", "Latitude"), crs = prj)
  }, error = function(e) {
    cat("Error when making species into points:", e$message, "\n")
  })

  if (verbose == T) {
    if(!any(is.na(sp_points))) cat("No values are NA:", green("TRUE") , "\n") else cat("Some values are NA:", red("FALSE") , "\n")   
  }
  
  return(sp_points)
}