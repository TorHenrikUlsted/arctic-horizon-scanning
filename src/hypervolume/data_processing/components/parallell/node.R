node_processing <- function(j, spec.list, proj.incl.t, method, accuracy, hv.projection, cores.max.high, min.disk.space, hv.dir, show.plot, verbose) {
  current_disk_space <- get_disk_space("/export", units = "GB")
  
  # Make dir objects
  log_dir <- paste0(hv.dir, "/logs")
  init_log <- paste0(log_dir, "/init-log.txt")
  create_file_if(init_log)
  node <- paste0("Node-", j)
  
  try(init_log <- file(init_log, open = "at"))
  sink(init_log, type = "output")
  sink(init_log, type = "message")
  
  node_timer <- start_timer(node)
  
  print(cores.max.high)
  
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
      break
    }
  }
  while_msg <- FALSE
  
  lock_node_it <- lock(lock_node_dir, lock.n = 1)
  
  print(lock_node_it)
  
  if (verbose) cat("Locked file. \n")
  
  node_con <- file(node_it, open = "a")
  identifier <- paste0("node", j)
  writeLines(identifier, node_it)
  close(node_con)
  
  unlock(lock_node_it)
  
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

  tryCatch({
    acq_data <- data_acquisition(show.plot, method, verbose = verbose, iteration = j, warn, err)
    
    invisible(gc())
    
    ana_data <- data_analysis(biovars_world = acq_data[[1]], biovars_region = acq_data[[2]], biovars_floreg = acq_data[[3]], method, show.plot, verbose = verbose, iteration = j, warn, err)
    
    invisible(gc())
    
    processed_data <- data_processing(spec, biovars_world = ana_data[[1]], spec.name, method, verbose = verbose, iteration = j, warn, err)
    invisible(gc())
    
    cat("Checking for locks... \n")
    
    while (TRUE) {
      if (is.locked(lock_hv_dir, lock.n = cores.max.high)) {
        if (!while_msg) {
          cat("The node is stuck in hypervolume analysis traffic... \n")
          while_msg <- TRUE
        }
        Sys.sleep(1)
      } else {
        Sys.sleep(runif(1, 0, 1))  # Add a random delay between 0 and 1 second
        break
      }
    }
    while_msg <- FALSE
    
    cat("After while loop. \n")
    
    lock_analysis <- lock(lock_hv_dir, lock.n = cores.max.high)
    
    analyzed_hv <- hv_analysis(processed_data, biovars_region = ana_data[[2]], region_hv = ana_data[[4]], method, spec.name, proj.incl.t, accuracy, hv.projection, verbose = verbose, iteration = j, proj_dir, warn, err)
    
    unlock(lock_analysis)
    
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
    
    fwrite(final_res, paste0(stats_dir, "/", method, "-stats.csv"), append = T, bom = T)
    
    while (TRUE) {
      if (!while_msg) {
        cat("The node is stuck in node-iterations traffic... \n")
        while_msg <- TRUE
      }
      if (is.locked(lock_node_dir, lock.n = 1)) {
        Sys.sleep(1)
      } else {
        Sys.sleep(runif(1, 0, 1))  # Add a random delay between 0 and 1 second
        break
      }
    }
    while_msg <- FALSE
    
    lock_node_it <- lock(lock_node_dir, lock.n = 1)
    
    node_con <- file(node_it, open = "r")
    lines <- readLines(node_it)
    close(node_con)
    lines <- lines[-which(grepl(identifier, lines))]
    node_con <- file(node_it, open = "w")
    writeLines(lines, node_it)
    close(node_con)
    
    unlock(lock_node_it)
    
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
