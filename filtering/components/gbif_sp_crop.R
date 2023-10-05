crop_species = function(sp_occ_df) {
  message("Starting cropping script")
  
  cat("No NAs in latitude: ", if(!any(is.na(sp_occ_df$decimalLatitude)) == T) {green(!any(is.na(sp_occ_df$decimalLatitude)))} else {red(any(is.na(sp_occ_df$decimalLatitude)))} ,"\n")
  cat("No NAs in longitude: ", if(!any(is.na(sp_occ_df$decimalLongitude)) == T) {green(!any(is.na(sp_occ_df$decimalLongitude)))} else {red(any(is.na(sp_occ_df$decimalLongitude)))} ,"\n")
  
  cat("Converting species to points \n")
  gbif_species_df_LongLat = cbind(sp_occ_df$decimalLongitude, sp_occ_df$decimalLatitude)
  sp_occ_points = vect(gbif_species_df_LongLat, sp_occ_df, type="points", crs = crs(cavm))
  
  cat("Cropping species to CAVM \n")
  sp_cavm_crop = crop(sp_occ, cavm)
  cat("Masking species to CAVM \n")
  sp_cavm_Mask = mask(sp_cavm_crop, cavm)
  plot(sp_cavm_Mask, add=T, col="red")
  
  cat("Making CAVM spatial points into dataframe and removing complete cases of NA \n")
  sp_occ_cavm = as.data.frame(sp_cavm_Mask)
  sp_occ_cavm = sp_occ_cavm[complete.cases(sp_occ_cavm), ]
  cat("All CAVM species are unique: ",  any(!unique(sp_occ_cavm$scientificName == "")) ,"\n")
  cat("No CAVM species are NA: ",  any(!is.na(sp_occ_cavm$scientificName == "")) ,"\n")
  
  sp_cavm = sp_occ_cavm[ , "scientificName", drop = F]
  cat(yellow("writing sp_cavm to csv: \n","./outputs/filtering/gbif_retrieval_process/gbif_cavm_species.csv \n"))
  write.csv(sp_cavm, "./outputs/filtering/gbif_retrieval_process/gbif_cavm_species.csv", row.names = F, fileEncoding = "UTF-8")
  
  cat("Cropping species to the Boreal region \n")
  sp_boreal_crop = crop(sp_occ, wwfecoRegions)
  cat("Masking species to the Boreal region \n")
  sp_boreal_Mask = mask(sp_boreal_crop, wwfecoRegions)
  plot(sp_boreal_Mask, add=T, col="green")
  
  cat("Making Boreal region spatial points into dataframe and removing complete cases of NA \n")
  sp_occ_boreal = as.data.frame(sp_boreal_Mask)
  sp_occ_boreal = sp_occ_boreal[complete.cases(sp_occ_boreal), ]
  cat("All CAVM species are unique: ",  any(!unique(sp_occ_boreal$scientificName == "")) ,"\n")
  cat("No CAVM species are NA: ",  any(!is.na(sp_occ_boreal$scientificName == "")) ,"\n")
  
  sp_boreal = sp_occ_boreal[ , "scientificName", drop = F]
  cat(yellow("writing sp_boreal to csv: \n","./outputs/filtering/gbif_retrieval_process/gbif_boreal_species.csv \n"))
  write.csv(sp_boreal, "./outputs/filtering/gbif_retrieval_process/gbif_boreal_species.csv", row.names = F, fileEncoding = "UTF-8")
  
  return(
    list(
      sp_cavm,
      sp_boreal
    )
  )
}