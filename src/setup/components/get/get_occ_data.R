get_occ_data <- function(species_w_keys, file.name, region = NULL, coord.uncertainty, download.key = NULL, download.doi = NULL) {
  vebcat("Initiating GBIF occurrence record download.", color = "funInit")

  download_path <- dirname(file.name)
  create_dir_if(download_path)

  species_keys <- species_w_keys$usageKey
  
  catn("Number of species keys:", highcat(length(species_keys)))
  vebprint(str(species_keys), text = "Species keys str:")
  
  out <- NULL

  if (is.null(download.key) & is.null(download.doi)) {
    catn(colcat("No", color = "indicator"), "download Key nor doi found.")
    
    if (is.null(region)) {
      catn("Region is", colcat("not", color = "indicator"), "being applied.")
    }  else {
      catn("Region", colcat("is", color = "indicator"), "being applied.")
    } 
    
    if (is.null(coord.uncertainty)) {
      catn("Coordinate Uncertainty is:", colcat("not", color = "indicator"), "being applied.")
    } else {
      catn("Coordinate Uncertainty is:", colcat(coord.uncertainty, color = "indicator"))
    }
    
    catn()

    u_input <- readline("Queue occurence download? Press [ENTER] to continue or [ESC] to exit.")

    if (u_input == "") {
      catn("Creating occucrence queue.")

      predicates <- list(
        pred_in("taxonKey", species_keys),
        pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN")),
        pred("hasCoordinate", TRUE),
        pred("hasGeospatialIssue", FALSE),
        pred("occurrenceStatus", "PRESENT"),
        pred_gte("year", 1970),
      )
      
      if (!is.null(coord.uncertainty)) {
        predicates <- c(predicates, list(pred_lte("coordinateUncertaintyInMeters", coord.uncertainty)))
      }
      
      if (!is.null(region)) {
        predicates <- c(predicates, list(pred("geometry", region)))
      }

      out <- do.call(occ_download, c(predicates, list(format = "SIMPLE_CSV")))

      tryCatch(
        {
          catn("Queueing download.")
          ## check status
          occ_download_wait(out)
          ## get the download Data and import to create dataframe
          gbif_occ_file <- occ_download_get(out, path = download_path, overwrite = T)
          
          download.key <- gbif_occ_file$keyname

          vebcat("GBIF occurrences Successfully downloaded.", color = "funSuccess")
        },
        error = function(e) {
          vebcat("An error has occurred: ", e$message, color = "nonFatalError")
        }
      )
    } else {
      stop("Process cancelled by user. \n")
    }
  } else {
    if (!is.null(download.key)) catn("Downlad key found.")
    if (!is.null(download.doi)) catn("Doi found.")
  }

  tryCatch(
    {
      if (!file.exists(paste0(file.name, ".csv"))) {
        if (file.exists(paste0(download_path, "/", download.key, ".zip"))) {
          catn("ZIP file named", paste0(download_path, "/", download.key, ".zip"), "found.")
        } else if (!file.exists(paste0(file.name, ".zip"))) {
          message("Trying to install by download key... ")
          catn("Using download key:", download.key)
          out <- occ_download_get(download.key, path = download_path, overwrite = TRUE)
        } else {
          catn("ZIP file named", paste0(download_path, "/", download.key, ".zip"), "found.")
        }

        catn("Unzipping GBIF file.")

        unzip(paste0(download_path, "/", download.key, ".zip"), exdir = download_path)

        csv <- paste0(download_path, "/", download.key, ".csv")
        
        catn("Renaming CSV file.")

        file.rename(from = csv, to = paste0(file.name, ".csv"))
        
        csv <- paste0(file.name, ".csv")

        catn("Reading GBIF CSV.")
        gbif_occ_df <- fread(csv)
        
        vebprint(head(gbif_occ_df, 3), text = "Sample of data table:")

        gbif_occ_df <- set_df_utf8(gbif_occ_df)

        fwrite(gbif_occ_df, paste0(file.name, ".csv"), bom = T)

        vebcat("GBIF occurrences Successfully downloaded.", color = "funSuccess")

        return(gbif_occ_df)
      } else {
        catn("GBIF Occurrence data found from", highcat(paste0(file.name, ".csv")))
        
        size <- file.size(paste0(file.name, ".csv")) / 1024^3 
        size <- round(size, digits = 2)
        catn("File size:", highcat(size), "GB.")
        
        if (size <= 5) {
          catn("File size smaller than 5 GB, reading file...")
          gbif_occ_df <- fread(paste0(file.name, ".csv"))
          return(gbif_occ_df)
        } else {
          vebcat("File size very big, using file chunking method.", color = "indicator")
        }
      }
    },
    error = function(e) { # 3
      vebcat("An error has occurred: ", e$message, color = "nonFatalError")

      message("Taking you to the online site... ")
      browseURL(download.doi)
    }
  )
}