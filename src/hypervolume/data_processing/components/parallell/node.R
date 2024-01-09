node_processing <- function(j, sp_list, proj.incl.t, method, accuracy, project, show.plot, verbose, min.disk.space) {
  current_disk_space <- get_disk_space("/home", units = "GB")
  
  warn <- function(w, warn_txt) {
    warn_msg <- conditionMessage(w)
    con <- file(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), open = "a")
    writeLines(paste(warn_txt, j, ":", warn_msg), con)
    close(con)
    invokeRestart(findRestart("muffleWarning"))
  }

  err <- function(e, err_txt) {
    err_msg <- conditionMessage(e)
    con <- file(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), open = "a")
    writeLines(paste(err_txt, j, ":", err_msg), con)
    close(con)
    stop(e)
  }

  # Check if the current disk space is less than the minimum required
  if (current_disk_space <= min.disk.space) stop("Insufficient disk space. Stopping processing.")
  
  print(sp_list[[j]])
  
  spec <- fread(sp_list[[j]], select = c(scientificName, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, coordinatePrecision, countryCode, stateProvince, year))
  
  print(names(spec))
  print(spec)

  spec.name <- names(spec)

  node_hv_dir <- paste0("./outputs/hypervolume/logs/nodes")
  create_dir_if(node_hv_dir)
  log_file <- paste0(node_hv_dir, "/", j, "-", gsub(" ", "-", spec.name), "-", method, "-log.txt")
  
  if (!file.exists(log_file)) {
    file.create(log_file)
  } else {
    file.remove(log_file)
    file.create(log_file)
  }

  try(log_file <- file(log_file, open = "at"))
  sink(log_file, type = "output")
  sink(log_file, type = "message")

  node_timer <- start_timer(paste0("Node-", j))

  cat("Run iteration", cc$lightSteelBlue(j), "\n")
  cat("Using species:", cc$lightSteelBlue(spec.name), "\n")

  tryCatch({
    acq_data <- data_acquisition(show.plot, method, verbose = verbose, iteration = j, warn, err)
    
    invisible(gc())
    
    ana_data <- data_analysis(biovars_world = acq_data[[1]], biovars_region = acq_data[[2]], biovars_floreg = acq_data[[3]], method, show.plot, verbose = verbose, iteration = j, warn, err)
    
    invisible(gc())
    
    processed_data <- data_processing(sp_list[[j]], biovars_world = ana_data[[1]], spec.name, method, projection, verbose = verbose, iteration = j, warn, err)
    invisible(gc())
    
    analyzed_hv <- hv_analysis(processed_data, biovars_region = ana_data[[2]], region_hv = ana_data[[4]], method, spec.name, proj.incl.t, accuracy = accuracy, project = project, verbose = verbose, iteration = j, warn, err)
    
    final_res <- data.frame(
      species = spec.name,
      observations = analyzed_hv[[1]],
      dimensions = analyzed_hv[[2]],
      samplesPerPoint = analyzed_hv[[3]],
      randomPoints = analyzed_hv[[4]],
      excluded = analyzed_hv[[5]],
      jaccard = analyzed_hv[[6]][[1]],
      sørensen = analyzed_hv[[6]][[2]],
      fracVolumeSpecies = analyzed_hv[[6]][[3]],
      fracVolumeRegion = analyzed_hv[[6]][[4]],
      overlapRegion = 1 - analyzed_hv[[6]][[4]],
      includedOverlap = sum(analyzed_hv[[7]] == T) / length(analyzed_hv[[7]])
    )
    # If we want the number of true and false values we can do: round(observations * includedOverlap) = true values, and: observations - true values = false values
    
    fwrite(final_res, paste0("./outputs/hypervolume/stats/", method, "-stats.csv"), append = T, bom = T)
    
    end_timer(node_timer)
    
    sink(type = "message")
    sink(type = "output")
    close(log_file)
    
  }, error = function(e) {
    sink(type = "message")
    sink(type = "output")
    close(log_file)
  })

  rm(list = setdiff(ls(), "j"))

  invisible(gc())

  list(
    iteration = j
  )
}
