occ_retrieval <- function(species_w_codes, region = NULL, download_key = NULL, doi = NULL, download_path, file_name) {
  message("Initiating GBIF Occurrence retrieval")
  species_codes <- species_w_codes$speciesKey
  sp_na <- any(is.na(species_codes))

  cat("Some codes are NA: ", if (sp_na == T) {
    red(sp_na)
    species_codes <- species_codes[!is.na(species_codes)]
  } else {
    green(sp_na)
  }, "\n Removing NAs... \n")
  cat("Some codes are blank: ", if (any(species_codes) == "") {
    cat(red("TRUE"), " \n Removing blanks... \n")
    species_codes <- species_codes[species_codes != ""]
  } else {
    cat(green("FALSE"))
  }, "\n")

  cat("Species codes str: ", str(species_codes), "\n")
  
  u_input <- readline("Queue occurence? Press [ENTER] to continue or [ESC] to exit.")

  if (u_input == "") {
    cat("Creating occucrence queue \n")
    
    out <- occ_download(
      pred_in("taxonKey", species_codes),
      pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN")),
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      if (!is.null(region)) pred("geometry", region),
      pred("occurrenceStatus", "PRESENT"),
      pred("year", 1970),
      pred_gte("coordinateUncertaintyInMeters", 5000),
      format = "SIMPLE_CSV"
    )
  } else {
    stop("Process cancelled by user. \n")
  }

  cat("Queueing download \n")

  tryCatch(
    {
      ## check status
      occ_download_wait(out)
      ## get the download Data and import to create dataframe
      gbif_occ_file <- occ_download_get(out, path = download_path, overwrite = T)

      cat("Renaming file to:", download_path, "\n")
      file.rename(from = gbif_occ_file, to = paste0(download_path, "/", file_name))
      gbif_occ_df <- occ_download_import(paste0(download_path, "/", file_name))
    },
    error = function(e) {
      cat(red("An error has occurred: ", e$message, "\n"))

      tryCatch(
        {
          message("Trying to install by download key... ")

          out <- occ_download_get(download_key, path = download_path, overwrite = TRUE)
          file.rename(from = out, to = file_name)
          gbif_occ_df <- occ_download_import(file_name)
        },
        error = function(e) {
          cat(red("An error has occurred: ", e$message, "\n"))

          message("Taking you to the online site... ")
          browseURL(doi)
        }
      )
    }
  )

  cat("Retrieved GBIF occurrences \n")

  return(gbif_occ_df)
}
