occ_retrieval = function(species_w_codes) {
  message("Initiating GBIF Occurrence retrieval")
  species_codes = species_w_codes$speciesKey
  sp_na = any(is.na(species_codes))
  
  cat("Some codes are NA: ", if(sp_na == T) {
    red(sp_na)
    species_codes = species_codes[!is.na(species_codes)]
  } else {green(sp_na)}, "\n Removing NAs... \n")
  cat("Some codes are blank: ", if (any(species_codes) == "") {
    cat(red("TRUE"), " \n Removing blanks... \n") 
    species_codes = species_codes[species_codes != ""]
  } else { cat(green("FALSE")) }, "\n")
  
  cat("Species codes str: ", str(species_codes), "\n")
  cat("Creating occucrence queue \n")
  
  out = occ_download(
    pred_in("taxonKey", species_codes),
    pred_in("basisOfRecord", c("HUMAN_OBSERVATION","PRESERVED_SPECIMEN")),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("geometry", combined_WKT),
    pred("occurrenceStatus","PRESENT"),
    format = "SIMPLE_CSV"
  )
  
  cat("Queueing download \n")
  
  tryCatch({
    ## check status
    occ_download_wait(out)
    ## get the download Data and import to create dataframe
    gbif_occ_file = occ_download_get(out, path = "resources", overwrite = T)
    
    cat("Renaming file to: resources/gbif_occ/gbif_occ_zip.zip \n")
    file.rename(from = gbif_occ_file, to = "resources/gbif_occ/gbif_occ_zip.zip")
    gbif_occ_df = occ_download_import("resources/gbif_occ/gbif_occ_zip.zip")
  }, error = function(e) {
    cat(red("An error has occurred: ", e$message, "\n"))
    
    tryCatch({
      message("Trying to install by DOI... ")
      
      doi = "https://doi.org/10.15468/dl.3dub7a"
      key = "0003202-230918134249559"
      new_path = "resources/gbif_occ/gbif_occ_zip.zip"
      out = occ_download_get(key, path = "resources", overwrite = TRUE)
      file.rename(from = out, to = new_path)
      gbif_occ_file = occ_download_import(out)
    }, error = function(e) {
      cat(red("An error has occurred: ", e$message, "\n"))
      
      message("Taking you to the online site... ")
      browseURL(doi)
    })
  })
  
  cat("Retrieved GBIF occurrences \n")
  
  return(gbif_occ_df)
}