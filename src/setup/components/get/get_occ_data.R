get_occ_data <- function(species_w_codes, download_path, file_name, region = NULL, download_key = NULL, doi = NULL) {
  cat(blue("Initiating GBIF occurrence record download. \n"))

  species_codes <- species_w_codes$speciesKey

  sp_na <- any(is.na(species_codes))

  if (sp_na == T) {
    cat(red("NA keys found, removing... \n"))
    species_codes <- species_codes[!is.na(species_codes)]
  } else if (any(species_codes) == "") {
    cat("Blank keys found, Removing... \n")
    species_codes <- species_codes[species_codes != ""]
  } else {
    cat(green("No blank keys, nor NAs found. \n"))
  }

  cat("Number of species:", length(species_codes), "\n")
  cat("Species codes str: \n")
  str(species_codes)

  out <- NULL

  if (is.null(download_key) & is.null(doi)) {
    cat("No download Key or doi found \n")

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

      tryCatch( # 1
        {
          if (!dir.exists(download_path)) dir.create(download_path, recursive = T)

          cat("Queueing download \n")

          ## check status
          occ_download_wait(out)
          ## get the download Data and import to create dataframe
          gbif_occ_file <- occ_download_get(out, path = download_path, overwrite = T)

          cat(cc$lightGreen("GBIF occurrences Successfully downloaded. \n"))

          return(gbif_occ_df)
        },
        error = function(e) {
          cat("Download failed with error. \n")
          cat(e, "\n")
          next
        }
      )
    } else {
      stop("Process cancelled by user. \n")
    }
  } else {
    if (!is.null(download_key)) cat("Downlad key found. \n")
    if (!is.null(doi)) cat("DOI found. \n")
  }

  if (!dir.exists(download_path)) dir.create(download_path, recursive = T)

  cat(red("An error has occurred: ", e$message, "\n"))

  tryCatch( # 2
    {
      message("Trying to install by download key... ")


      cat("Using download key:", download_key, "\n")

      if (!file.exists(paste0(file_name, ".csv"))) {
        out <- occ_download_get(download_key, path = download_path, overwrite = TRUE)

        cat("Renaming file from:", out, "to:", paste0(file_name, ".zip"), "\n")

        if (file.exists(paste0(file_name, ".zip"))) {
          cat("File already exists. \n")
        } else if (!file.rename(from = out, to = paste0(file_name, ".zip"))) {
          cat("Failed to rename the file. Please check the file paths and permissions.\n")
        } else {
          cat("File renamed successfully.\n")
        }

        cat("Unzipping GBIF file. \n")

        unzip(paste0(file_name, ".zip"), exdir = download_path)

        csv <- paste0(download_path, "/", download_key, ".csv")

        file.rename(from = csv, to = paste0(file_name, ".csv"))

        csv <- paste0(file_name, ".csv")

        cat("Reading GBIF CSV. \n")
        gbif_occ_df <- fread(csv, sep = "\t")

        cat("Sample of data table: \n")
        print(head(gbif_occ_df, 3))

        gbif_occ_df <- set_df_utf8(gbif_occ_df)

        fwrite(gbif_occ_df, paste0(file_name, ".csv"), bom = T)

        cat(cc$lightGreen("GBIF occurrences Successfully downloaded. \n"))

        return(gbif_occ_df)
      }
    },
    error = function(e) { # 3
      cat(red("An error has occurred: ", e$message, "\n"))

      message("Taking you to the online site... ")
      browseURL(doi)
    }
  )
}
