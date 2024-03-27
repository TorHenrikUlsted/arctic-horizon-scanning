get_sp_list <- function(taxon, region, region.name, file.name, download.key, download.doi) {
  vebcat("Initiating download of GBIF species list based on region.", color = "funInit")

  download_path <- dirname(file.name)
  create_dir_if(download_path)

  catn("Getting species within the", highcat(taxon), "taxon.")

  taxon_key <- name_backbone(taxon)$usageKey

  if(!is.null(region.name)) {
    catn("Applying the", highcat(region.name), "region.")
  } else {
    catn("Region will NOT be applied.")
  }

  out <- NULL

  if (is.null(download.key) & is.null(download.doi)) {
    catn("No download Key or download.doi found")

    u_input <- readline("Queue occurence download? Press [ENTER] to continue or [ESC] to exit.")

    if (u_input == "") {
      catn("Creating occucrence queue.")

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
          catn("Queueing download.")
          ## check status
          occ_download_wait(out)
          ## get the download Data and import to create dataframe
          gbif_sp_list <- occ_download_get(out, path = download_path, overwrite = T)
        },
        error = function(e) {
          vebcat("Error when trying to get download.", color = "fatalError")
          stop(e)
        }
      )
    } else {
      stop("Process cancelled by user.")
    }
  } else {
    if (!is.null(download.key)) catn("Downlad key found.")
    if (!is.null(download.doi)) catn("download.doi found.")
    
    file_key <- paste0(download_path, "/", download.key, ".zip")

    tryCatch(
      {
        if (!file.exists(paste0(file.name, ".csv"))) {
          catn("CSV file not found.")
          if (!file.exists(file_key)) {
            catn("ZIP file not found.")
            catn("Trying to install using download key", highcat(download.key))
            out <- occ_download_get(download.key, path = download_path, overwrite = TRUE)
          }
        } else {
          catn("CSV file found.")
          vebcat("GBIF species list downloaded successfully.", funSuccess)
          return(gbif_sp_list = fread(paste0(file.name, ".csv")))
        }
        
      },
      error = function(e) {
        catn("Error when trying to install using download key.")
        stop(e)
      }
    )
  }

  tryCatch(
    {
      if (!file.exists(paste0(file.name, ".csv"))) {
        
        catn("Unzipping GBIF file.")

        unzip(file_key, exdir = download_path)

        csv <- paste0(download_path, "/", download.key, ".csv")

        file.rename(from = csv, to = paste0(file.name, ".csv"))

        csv <- paste0(file.name, ".csv")

        catn("Reading GBIF CSV.")
        gbif_sp_list <- fread(csv, sep = "\t")

        catn("Sample of data table:")
        print(head(gbif_sp_list, 3))

        gbif_sp_list <- set_df_utf8(gbif_sp_list)

        fwrite(gbif_sp_list, paste0(file.name, ".csv"), bom = T)

        vebcat("GBIF species list downloaded successfully.", color = "funSuccess")

        return(gbif_sp_list)
      }
    },
    error = function(e) { # 3
      vebcat("Error when unzipping GBIF ZIP file:", e$message, color = "nonFatalError")
      
      if (!file.exists(paste0(file.name, ".zip"))) {
        catn("GBIF ZIP file not found.")
        message("Taking you to the online site... ")
        Sys.sleep(3)
        browseURL(download.doi)
      }
    }
  )
}
