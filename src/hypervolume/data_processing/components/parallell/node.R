node_processing <- function(j, spec.list, proj.incl.t, method, accuracy, ndim, hv.projection, cores.max.high, min.disk.space, hv.dir, show.plot, verbose) {
  current_disk_space <- get_disk_space("/export", units = "GB")
  skip_to_end <- FALSE
  
  # Make dir objects
  log_dir <- paste0(hv.dir, "/logs")
  init_log <- paste0(log_dir, "/init-log.txt")
  create_file_if(init_log)
  node <- paste0("Node-", j)
  
  try(init_log <- file(init_log, open = "at")) # for error handling before the process even starts
  sink(init_log, type = "output")
  sink(init_log, type = "message")
  
  node_timer <- start_timer(node)
  
  cat("Cores.max.high:", cores.max.high, "\n")
  
  if (verbose) cat("Creating directories. \n")
  
  node_hv_dir <- paste0(log_dir, "/nodes")
  stats_dir <- paste0(hv.dir, "/stats")
  proj_dir <- paste0(hv.dir, "/projections")
  locks_dir <- paste0(hv.dir, "/locks")
  lock_node_dir <- paste0(locks_dir, "/node-iteration")
  lock_hv_dir <- paste0(locks_dir, "/hv_analysis")
  
  # Create dirs of they do not exist
  create_dir_if(node_hv_dir)
  create_dir_if(lock_node_dir)
  create_dir_if(lock_hv_dir)
  
  if (verbose) cat("Creating files. \n")
  
  # Make file objects
  warn_file <- paste0(log_dir, "/", method, "-warning.txt")
  err_file <- paste0(log_dir, "/", method, "-error.txt")
  node_it <- paste0(log_dir, "/node-iterations.txt")
  while_msg <- FALSE
  
  if (verbose) cat("Locking file. \n")
  
  while (TRUE) {
    if (!while_msg) {
      cat("The node is stuck in node-iterations traffic... \n")
      while_msg <- TRUE
    }
    if (is.locked(lock_node_dir, lock.n = 1)) {
      Sys.sleep(1)
    } else {
      Sys.sleep(runif(1, 0, 1))  # Add a random delay between 0 and 1 second
      lock_node_it <- lock(lock_node_dir, lock.n = 1)
      if (!is.null(lock_node_it) && file.exists(lock_node_it)) {
        break
      }
    }
  }
  while_msg <- FALSE
  
  if (verbose) cat("Checking for locked file. \n")
  
  node_con <- file(node_it, open = "at")
  identifier <- paste0("node", j)
  writeLines(identifier, node_con)
  close(node_con)
  
  unlock(lock_node_it)
  
  cat("spec.list[[j]]:\n")
  print(spec.list[[j]])
  
  # Check if the current disk space is less than the minimum required
  if (current_disk_space <= min.disk.space) stop("Insufficient disk space. Stopping processing.")
  
  spec_to_read <- spec.list[[j]]
  cat("spec_to_read:\n")
  print(spec_to_read)
  
  spec.name <- gsub(".csv", "", basename(spec_to_read))
  spec.name <- trimws(spec.name)
  cat("spec.name:\n", spec.name, "\n")
  
  
  spec <- fread(spec_to_read, select = c("species", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", "coordinatePrecision", "countryCode", "stateProvince", "year"))
  
  nobs <- nrow(spec)
  
  cat("spec:\n")
  print(head(spec, 2))

  cat("Init-log finished. \n")
  
  sink(type = "message")
  sink(type = "output")
  close(init_log)
  

  main_log <- paste0(node_hv_dir, "/", j, "-", gsub(" ", "-", spec.name), "-", method, "-log.txt")
  create_file_if(main_log)
  
  try(main_log <- file(main_log, open = "at"))
  sink(main_log, type = "output")
  sink(main_log, type = "message")

  cat("Run iteration", cc$lightSteelBlue(j), "\n")
  cat("Using species:", cc$lightSteelBlue(spec.name), "\n")
  cat("Species observations:", nobs, "\n")
  cat("Log observations:", log(nobs), "\n")
  cat("Expected dimensions:", ndim, "\n")
  cat("Identifier:", identifier, "\n\n")
  
  if (log(nobs) < ndim) {
    cat("Number of observations are less than the expected dimensions. \n")
    cat("Appending data to csv file. \n")
    
    final_res <- data.frame(
      species = spec.name,
      observations = nobs,
      dimensions = ndim,
      samplesPerPoint = 0,
      randomPoints = 0,
      excluded = T,
      jaccard = 0,
      sørensen = 0,
      fracVolumeSpecies = 0,
      fracVolumeRegion = 0,
      overlapRegion = 1,
      includedOverlap = 0
    )
    
    skip_to_end <- TRUE
  }

  tryCatch({
    if (verbose) cat("Creating warning and error functions. \n")
    # Create warning and error functions
    warn <- function(w, warn_txt) {
      warn_msg <- conditionMessage(w)
      warn_con <- file(warn_file, open = "a")
      writeLines(paste(warn_txt, j, ":", warn_msg), warn_con)
      close(warn_con)
      invokeRestart(findRestart("muffleWarning"))
    }
    
    err <- function(e, err_txt) {
      err_msg <- conditionMessage(e)
      err_con <- file(err_file, open = "a")
      writeLines(paste(err_txt, j, ":", err_msg), err_con)
      close(err_con)
      stop(e)
    }
    
  if (!skip_to_end) {
    cat("Running main sequence. \n")
    
    acq_data <- data_acquisition(show.plot, method, verbose = verbose, iteration = j, warn = warn, err = err)
    
    invisible(gc())
    
    ana_data <- data_analysis(biovars_world = acq_data[[1]], biovars_region = acq_data[[2]], biovars_floreg = acq_data[[3]], method, show.plot, verbose = verbose, iteration = j, warn = warn, err = err)
    
    invisible(gc())
    
    processed_data <- data_processing(spec, biovars_world = ana_data[[1]], spec.name, method, verbose = verbose, iteration = j, warn = warn, err = err)
    invisible(gc())
    
    
    analyzed_hv <- hv_analysis(processed_data, biovars_region = ana_data[[2]], region_hv = ana_data[[4]], method, spec.name, proj.incl.t, accuracy, hv.projection, verbose = verbose, iteration = j, proj_dir, lock_hv_dir, cores.max.high, warn = warn, err = err)
    
    cat("Appending data to csv file. \n")
    
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
  } else {
    cat("Skipping to end. \n")
  }
    
    fwrite(final_res, paste0(stats_dir, "/", method, "-stats.csv"), append = T, bom = T)
    
    cat("Appended data. \n")
    print(final_res)
    
    while (TRUE) {
      if (!while_msg) {
        cat("The node has entered node-iterations queue... \n")
        while_msg <- TRUE
      }
      if (is.locked(lock_node_dir, lock.n = 1)) {
        Sys.sleep(1)
      } else {
        Sys.sleep(runif(1, 0, 1))  # Add a random delay between 0 and 1 second
        lock_node_it <- lock(lock_node_dir, lock.n = 1)
        if (!is.null(lock_node_it) && file.exists(lock_node_it)) {
          break
        }
      }
    }
    while_msg <- FALSE
    
    node_con <- file(node_it, open = "r")
    lines <- readLines(node_con)
    close(node_con)
    if (verbose) cat("Removing from node-iterations:", lines[which(grepl(paste0("^", identifier, "$"), lines))], "\n")
    lines <- lines[-which(grepl(paste0("^", identifier, "$"), lines))] # Use regular expression to get exact match
    if (verbose) cat("Rewriting to node-iterations:", lines, "\n")
    node_con <- file(node_it, open = "w")
    writeLines(lines, node_con)
    close(node_con)
    
    cat("Unlocking node iterations. \n")
    
    unlock(lock_node_it)
    
    cat("Node finished successfully. \n")
    
    end_timer(node_timer)
    
    sink(type = "message")
    sink(type = "output")
    close(main_log)
    
  }, error = function(e) {
    cat("Error occurred in node", j, ":", e$message, "\n")
    sink(type = "message")
    sink(type = "output")
    close(main_log)
  })

  rm(list = setdiff(ls(), "j"))

  invisible(gc())

  list(
    iteration = j
  )
}
