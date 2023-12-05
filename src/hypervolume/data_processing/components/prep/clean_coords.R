clean_coords <- function(df, verbose = T) {
  if (!is.data.frame(df)) {
    stop("Input must be a data.frame or data.table")
  }
  
  if (verbose == T) cat("Getting Long/Lat values. \n")
  
  df <- df %>% 
    distinct(decimalLongitude, decimalLatitude, species, .keep_all = TRUE)
  
  any(duplicated(paste(df$decimalLongitude, df$decimalLatitude)))
  
  
  cat("Cleaning species using coordinateCleaner. \n")
  create_dir_if("./outputs/data_processing/prep")
  # "flag" sets adds true/false to the corresponding tests
  df <- clean_coordinates(df, value = "clean")
  
  cat("Thining coordinates with spThin. \n")
  
  sp_thinned <- thin(
    loc.data = df,
    lat.col = "decimalLatitude",
    long.col = "decimalLongitude",
    spec.col = "species",
    write.files = T,
    max.files = 1,
    out.dir = "./outputs/data_processing/prep",
    out.base = "sp",
    thin.par = 0.1,
    reps = 10,
  )
  
  return(sp_thinned)
}