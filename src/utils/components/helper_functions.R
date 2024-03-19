##########################
#       Data.table       #
##########################

merge_and_sum <- function(dt1, dt2, sumCol, by, all = TRUE) {
  merged_dt <- merge(dt1, dt2, by = by, all = all)
  merged_dt[[sumCol]] <- rowSums(merged_dt[, c(paste0(sumCol,".x"), paste0(sumCol, ".y"))], na.rm = TRUE)
  
  return(merged_dt)
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
