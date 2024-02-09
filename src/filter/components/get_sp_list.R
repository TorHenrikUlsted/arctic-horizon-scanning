get_sp_list <- function(taxon, region, region.name, file.name, download.key, download.doi) {
  cat(blue("Initiating download of GBIF species list based on region. \n"))

  download_path <- dirname(file.name)
  create_dir_if(download_path)

  cat("Getting species within the", cc$lightSteelBlue(taxon), "taxon. \n")

  taxon_key <- name_backbone(taxon)$usageKey

  cat("Applying the", cc$lightSteelBlue(region.name), "region. \n")

  out <- NULL

  if (is.null(download.key) & is.null(download.doi)) {
    cat("No download Key or download.doi found \n")

    u_input <- readline("Queue occurence download? Press [ENTER] to continue or [ESC] to exit.")

    if (u_input == "") {
      cat("Creating occucrence queue \n")

      predicates <- list(
        pred_in("taxonKey", taxon_key),
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

      out <- do.call(occ_download, c(predicates, list(format = "SPECIES_LIST")))

      tryCatch(
        {
          cat("Queueing download \n")
          ## check status
          occ_download_wait(out)
          ## get the download Data and import to create dataframe
          gbif_sp_list <- occ_download_get(out, path = download_path, overwrite = T)
        },
        error = function(e) {
          cat(red("Error when trying to get download. \n"))
          stop(e)
        }
      )
    } else {
      stop("Process cancelled by user. \n")
    }
  } else {
    if (!is.null(download.key)) cat("Downlad key found. \n")
    if (!is.null(download.doi)) cat("download.doi found. \n")
    
    file_key <- paste0(download_path, "/", download.key, ".zip")

    tryCatch(
      {
        if (!file.exists(paste0(file.name, ".csv"))) {
          cat("CSV file not found. \n")
          if (!file.exists(file_key)) {
            cat("ZIP file not found. \n")
            cat("Trying to install using download key", cc$lightSteelBlue(download.key), "\n")
            out <- occ_download_get(download.key, path = download_path, overwrite = TRUE)
          }
        } else {
          cat("CSV file found. \n")
          cat(cc$lightGreen("GBIF species list downloaded successfully. \n"))
          return(gbif_sp_list = fread(paste0(file.name, ".csv")))
        }
        
      },
      error = function(e) {
        cat("Error when trying to install using download key. \n")
        stop(e)
      }
    )
  }

  tryCatch(
    {
      if (!file.exists(paste0(file.name, ".csv"))) {
        
        cat("Unzipping GBIF file. \n")

        unzip(file_key, exdir = download_path)

        csv <- paste0(download_path, "/", download.key, ".csv")

        file.rename(from = csv, to = paste0(file.name, ".csv"))

        csv <- paste0(file.name, ".csv")

        cat("Reading GBIF CSV. \n")
        gbif_sp_list <- fread(csv, sep = "\t")

        cat("Sample of data table: \n")
        print(head(gbif_sp_list, 3))

        gbif_sp_list <- set_df_utf8(gbif_sp_list)

        fwrite(gbif_sp_list, paste0(file.name, ".csv"), bom = T)

        cat(cc$lightGreen("GBIF species list downloaded successfully. \n"))

        return(gbif_sp_list)
      }
    },
    error = function(e) { # 3
      cat(red("Error when unzipping GBIF ZIP file:", e$message, "\n"))
      
      if (!file.exists(paste0(file.name, ".zip"))) {
        cat("GBIF ZIP file not found. \n")
        message("Taking you to the online site... ")
        Sys.sleep(3)
        browseURL(download.doi)
      }
    }
  )
}
