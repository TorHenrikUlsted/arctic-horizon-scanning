##########################
#       Data.table       #
##########################

merge_and_sum <- function(dt1, dt2, sumCol, by, all = TRUE) {
  merged_dt <- merge(dt1, dt2, by = by, all = all)
  merged_dt[[sumCol]] <- rowSums(merged_dt[, c(paste0(sumCol,".x"), paste0(sumCol, ".y"))], na.rm = TRUE)
  
  return(merged_dt)
}

find_min_data_col <- function(file_path, verbose = TRUE) {
  catn("Finding column with the least memory allocation needed.")
  # Read the first 10 rows of the file
  dt <- fread(file_path, nrows = 10)
  
  # Calculate the total character length of each column
  column_lengths <- sapply(dt, function(x) sum(nchar(as.character(x))))
  
  # Find the column name with the least character length
  least_data_column <- names(column_lengths)[which.min(column_lengths)]
  
  catn("Column with the least data:", highcat(least_data_column))
  
  return(least_data_column)
}

extract_name_after_prefix = function(x, prefixes = c("var", "subsp.", "ssp.", "f.")) {
  
  # Create a regular expression pattern to match the specified prefixes
  prefix_pattern = paste0(prefixes, collapse = "|")
  
  # Use str_extract to extract the name after the specified prefixes
  name = str_extract(x, paste0("(?<=", prefix_pattern, ")\\s*\\S+"))
  
  return(name)
} 

log_duplicates <- function(df, column, process, folder_name, file_name) {
  
  df <- as.data.frame(df)
  
  dup_rows <- which(duplicated(df[[column]]))
  dup_sp <- df[dup_rows, column]
  
  log_out_df <- data.frame(scientificName = dup_sp, rowNumber = dup_rows)
  
  dups_n <- nrow(log_out_df)
  
  log_path <- file.path("./outputs", process, "logs", folder_name)
  create_dir_if(log_path)
  
  output_log_path <- file.path(log_path, file_name)
  
  fwrite(log_out_df, file = output_log_path, bom = T)
  
  catn("Number of duplicates:", cc$aquamarine(dups_n))
  catn("More info found in: ", yellow(output_log_path))
}

set_df_utf8 <- function(df) {
  
  for (name in names(df)[sapply(df, is.character)]) {
    df[[name]] <- enc2utf8(df[[name]])
  }
  
  return(df)
}

##########################
#         Terra          #
##########################

extract_ext_to_dt <- function(raster, value = "value", cells = TRUE) {
  rast_extr <- terra::extract(raster, ext(raster), cells = cells)
  rast_dt <- as.data.table(rast_extr)
  names(rast_dt) <- c("cell", value)
  
  return(rast_dt)
}

convert_spatial_dt <- function(spatial, verbose = FALSE) {
  dt <- as.data.frame(spatial)
  dt <- as.data.table(dt)
  
  return(dt)
}

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

order_by_apg <- function(input, by, verbose = FALSE) {
  # Load apg4
  apg <- fread("./resources/taxon/apg4/apg4.txt")
  
  apg_rank <- apg[taxonRank == by]
  
  apg_vect <- apg_rank$scientificName
  
  non_apg <- setdiff(input, apg_vect)
  
  vebcat("non_apg:", veb = verbose)
  vebprint(non_apg, veb = verbose)
  
  input_apg <- input[!input %in% non_apg]

  vebcat("input without non_apgs:", veb = verbose)
  vebprint(input_apg, veb = verbose)
  
  vebcat("Ordering apgs", veb = verbose)
  if (is.vector(input)) {
    ordered_apg <- apg_vect[apg_vect %in% input_apg]
    
    # Add the others
    non_apg_order <- c("Lycopodiales", "Selaginellales", "Equisetales", "Osmundales", "Salviniales", "Cyatheales", "Polypodiales","Ephedrales","Pinales")
    
    ordered <- c(non_apg_order, ordered_apg)
    print(ordered)
    
  } else {
    setorderv(input, cols = apg_rank[[by]])
  }
  
  return(ordered)
}

find_peaks <- function(data, prominence = 0.1) {
  # Identify local maxima using diff
  peaks <- which(diff(data) > 0 & diff(c(data, 0)) < 0)
  
  # Filter based on prominence (difference from neighbors)
  filtered_peaks <- peaks[data[peaks] - data[peaks - 1] >= prominence & data[peaks] - data[peaks + 1] >= prominence]
  
  return(filtered_peaks)
}

##########################
#       Parallel         #
##########################
load_sp_rast <- function(spec.filename) {
  sp_name <- basename(dirname(spec.filename))
  
  sp_rast <- terra::rast(spec.filename)
  names(sp_rast) <- sp_name
  
  return(sp_rast)
}

##########################
#      File system       #
##########################

source_all <- function(dir) {
  # Get a list of all .R files in the directory and its subdirectories
  r_files <- list.files(path = dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  
  # Source each file
  lapply(r_files, source)
  
  cat(length(r_files), "scripts sourced from", dir, "and its subdirectories.\n")
}

