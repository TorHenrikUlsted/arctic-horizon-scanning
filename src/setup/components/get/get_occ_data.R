get_occ_data <- function(species_w_keys, file.name, region = NULL, download.key = NULL, download.doi = NULL) {
  cat(blue("Initiating GBIF occurrence record download. \n"))

  download_path <- dirname(file.name)
  create_dir_if(download_path)

  species_codes <- species_w_keys$speciesKey

  sp_na <- any(is.na(species_codes))

  if (sp_na == T) {
    cat(red("NA keys found, removing... \n"))
    species_codes <- species_codes[!is.na(species_codes)]
  } else if (any(species_codes) == "") {
    cat("Blank keys found, Removing... \n")
    species_codes <- species_codes[species_codes != ""]
  } else {
    cat(cc$lightGreen("No blank keys, nor NAs found. \n"))
  }

  cat("Number of species:", cc$lightSteelBlue(length(species_codes)), "\n")
  cat("Species codes str: \n")
  str(species_codes)

  out <- NULL

  if (is.null(download.key) & is.null(download.doi)) {
    cat(cc$aquamarine("No"), "download Key nor doi found \n")
    if (is.null(region)) cat("Region is", cc$aquamarine("not"), "being applied. \n") else cat("Region", cc$aquamarine("is"), "being applied. \n")

    u_input <- readline("Queue occurence download? Press [ENTER] to continue or [ESC] to exit.")

    if (u_input == "") {
      cat("Creating occucrence queue \n")

      predicates <- list(
        pred_in("taxonKey", species_codes),
        pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "PRESERVED_SPECIMEN")),
        pred("hasCoordinate", TRUE),
        pred("hasGeospatialIssue", FALSE),
        pred("occurrenceStatus", "PRESENT"),
        pred_gte("year", 1970),
        pred_lte("coordinateUncertaintyInMeters", 4700) # 2.5 min res = 4.63 km at equator, WorldClim Cell size in km = 111.32 * cos(60 * pi/180) * 2.5 / 60 ≈ 2.8 km in the boreal belt and 1.9 km in the arctic.
      )

      if (!is.null(region)) {
        predicates <- c(predicates, list(pred("geometry", region)))
      }

      out <- do.call(occ_download, c(predicates, list(format = "SIMPLE_CSV")))

      tryCatch(
        {
          cat("Queueing download \n")
          ## check status
          occ_download_wait(out)
          ## get the download Data and import to create dataframe
          gbif_occ_file <- occ_download_get(out, path = download_path, overwrite = T)

          cat(cc$lightGreen("GBIF occurrences Successfully downloaded. \n"))

          return(gbif_occ_file)
        },
        error = function(e) {
          cat(red("An error has occurred: ", e$message, "\n"))
        }
      )
    } else {
      stop("Process cancelled by user. \n")
    }
  } else {
    if (!is.null(download.key)) cat("Downlad key found. \n")
    if (!is.null(download.doi)) cat("Doi found. \n")
  }

  tryCatch(
    {
      if (!file.exists(paste0(file.name, ".csv"))) {
        if (file.exists(paste0(download_path, "/", download.key, ".zip"))) {
          cat("ZIP file named", paste0(download_path, "/", download.key, ".zip"), "found. \n")
          if (!file.rename(from = out, to = paste0(file.name, ".zip"))) {
            cat("Failed to rename the file. Please check the file paths and permissions.\n")
          } else {
            cat("File renamed successfully.\n")
          }
        } else if (!file.exists(paste0(file.name, ".zip"))) {
          message("Trying to install by download key... ")
          cat("Using download key:", download.key, "\n")
          out <- occ_download_get(download.key, path = download_path, overwrite = TRUE)
          if (!file.rename(from = out, to = paste0(file.name, ".zip"))) {
            cat("Failed to rename the file. Please check the file paths and permissions.\n")
          } else {
            cat("File renamed successfully.\n")
          }
        } else {
          cat("ZIP file named", paste0(file.name, ".zip"), "found. \n")
        }

        cat("Unzipping GBIF file. \n")

        unzip(paste0(file.name, ".zip"), exdir = download_path)

        csv <- paste0(download_path, "/", download.key, ".csv")
        
        cat("Renaming CSV file. \n")

        file.rename(from = csv, to = paste0(file.name, ".csv"))

        csv <- paste0(file.name, ".csv")

        cat("Reading GBIF CSV. \n")
        gbif_occ_df <- fread(csv, sep = "\t")

        cat("Sample of data table: \n")
        print(head(gbif_occ_df, 3))

        gbif_occ_df <- set_df_utf8(gbif_occ_df)

        fwrite(gbif_occ_df, paste0(file.name, ".csv"), bom = T)

        cat(cc$lightGreen("GBIF occurrences Successfully downloaded. \n"))

        return(gbif_occ_df)
      } else {
        cat("GBIF Occurrence data found from", cc$lightSteelBlue(paste0(file.name, ".csv")), ".\n")
        
        size <- file.size(paste0(file.name, ".csv")) / 1024^3 
        size <- round(size, digits = 2)
        cat("File size:", cc$lightSteelBlue(size), "GB. \n")
        
        if (size <= 10) {
          gbif_occ_df <- fread(paste0(file.name, ".csv"), sep = "\t")
          return(gbif_occ_df)
        } else {
          cat(cc$aquamarine("File size very big, load manually. \n"))
        }
      }
    },
    error = function(e) { # 3
      cat(red("An error has occurred: ", e$message, "\n"))

      message("Taking you to the online site... ")
      browseURL(download.doi)
    }
  )
}