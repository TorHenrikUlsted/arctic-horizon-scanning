node_processing <- function(j, spec.list, proj.incl.t, method, accuracy, hv.dims, hv.projection, cores.max.high, min.disk.space, hv.dir, show.plot, verbose) {
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

  catn("Cores.max.high:", cores.max.high)

  vebcat("Creating directories.", veb = verbose)

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

  vebcat("Creating files.", veb = verbose)

  # Make file objects
  warn_file <- paste0(log_dir, "/", method, "-warning.txt")
  err_file <- paste0(log_dir, "/", method, "-error.txt")
  node_it <- paste0(log_dir, "/node-iterations.txt")
  while_msg <- FALSE

  vebcat("Locking file.", veb = verbose)

  while (TRUE) {
    if (!while_msg) {
      vebcat("The node is stuck in node-iterations traffic...", veb = verbose)
      while_msg <- TRUE
    }
    if (is.locked(lock_node_dir, lock.n = 1)) {
      Sys.sleep(1)
    } else {
      Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
      lock_node_it <- lock(lock_node_dir, lock.n = 1)
      if (!is.null(lock_node_it) && file.exists(lock_node_it)) {
        break
      }
    }
  }
  while_msg <- FALSE

  vebcat("Checking for locked file.", veb = verbose)

  node_con <- file(node_it, open = "at")
  identifier <- paste0("node", j)
  writeLines(identifier, node_con)
  close(node_con)

  unlock(lock_node_it)

  vebprint(spec.list[[j]], text = "spec.list[[j]]:")

  # Check if the current disk space is less than the minimum required
  if (current_disk_space <= min.disk.space) stop("Insufficient disk space. Stopping processing.")

  spec_to_read <- spec.list[[j]]
  vebprint(spec_to_read, text = "spec_to_read")
  
  spec.name <- gsub(".csv", "", basename(spec_to_read))
  spec.name <- trimws(spec.name)
  catn("spec.name:\n", spec.name)

  spec <- fread(spec_to_read, select = c("species", "decimalLongitude", "decimalLatitude", "coordinateUncertaintyInMeters", "coordinatePrecision", "countryCode", "stateProvince", "year"))

  nobs <- nrow(spec)

  vebprint(head(spec, 2), text = "spec:")

  catn("Init-log finished.")

  sink(type = "message")
  sink(type = "output")
  close(init_log)

  main_log <- paste0(node_hv_dir, "/", j, "-", gsub(" ", "-", spec.name), "-", method, "-log.txt")
  create_file_if(main_log)

  try(main_log <- file(main_log, open = "at"))
  sink(main_log, type = "output")
  sink(main_log, type = "message")

  catn("Run iteration", highcat(j))
  catn("Using species:", highcat(spec.name))
  catn("Species observations:", nobs)
  catn("Log observations:", log(nobs))
  catn("Expected dimensions:", length(hv.dims))
  catn("Identifier:", identifier, "\n")

  if (log(nobs) < length(hv.dims)) {
    catn("Number of observations are less than the expected dimensions.")
    skip_to_end <- TRUE
  }

  tryCatch(
    {
      vebcat("Creating warning and error functions.", veb = verbose)
      # Create warning and error functions
      warn <- function(w, warn_txt) {
        warn_msg <- conditionMessage(w)
        warn_con <- file(warn_file, open = "a")
        writeLines(paste(warn_txt, "in iteration", j, ":", warn_msg), warn_con)
        close(warn_con)
        invokeRestart(findRestart("muffleWarning"))
      }

      err <- function(e, err_txt) {
        err_msg <- conditionMessage(e)
        err_con <- file(err_file, open = "a")
        writeLines(paste(err_txt, "in iteration", j, ":", err_msg), err_con)
        close(err_con)
        stop(e)
      }

      while (!skip_to_end) {
        catn("Running main sequence.")

        acq_data <- data_acquisition(show.plot, method, verbose = verbose, iteration = j, warn = warn, err = err)

        invisible(gc())

        ana_data <- data_analysis(biovars_world = acq_data[[1]], biovars_region = acq_data[[2]], biovars_floreg = acq_data[[3]], hv.dims, method, show.plot, verbose = verbose, iteration = j, warn = warn, err = err)

        invisible(gc())

        processed_data <- data_processing(spec, biovars_world = ana_data[[1]], spec.name, method, verbose = verbose, iteration = j, warn = warn, err = err)
        
        invisible(gc())

        if (is.list(processed_data) && processed_data$excluded == TRUE) {
          skip_to_end <- TRUE
          break
        }

        analyzed_hv <- hv_analysis(processed_data, biovars_region = ana_data[[2]], region_hv = ana_data[[4]], method, spec.name, proj.incl.t, accuracy, hv.projection, verbose = verbose, iteration = j, proj_dir, lock_hv_dir, cores.max.high, warn = warn, err = err)

        catn("Appending data to csv file.")

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
        
        break
      }

      if (skip_to_end) {
        catn("Skipping to end.")
        catn("Appending data to csv file.")

        final_res <- data.frame(
          species = spec.name,
          observations = nobs,
          dimensions = length(hv.dims),
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
      }

      fwrite(final_res, paste0(stats_dir, "/", method, "-stats.csv"), append = T, bom = T)

      vebprint(final_res, text = "Appended data:")

      while (TRUE) {
        if (!while_msg) {
          catn("The node has entered node-iterations queue...")
          while_msg <- TRUE
        }
        if (is.locked(lock_node_dir, lock.n = 1)) {
          Sys.sleep(1)
        } else {
          Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
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
      vebcat("Removing from node-iterations:", lines[which(grepl(paste0("^", identifier, "$"), lines))], veb = verbose)
      lines <- lines[-which(grepl(paste0("^", identifier, "$"), lines))] # Use regular expression to get exact match
      vebcat("Rewriting to node-iterations:", lines, veb = verbose)
      node_con <- file(node_it, open = "w")
      writeLines(lines, node_con)
      close(node_con)

      catn("Unlocking node iterations.")

      unlock(lock_node_it)

      catn("Node finished successfully.")

      end_timer(node_timer)

      sink(type = "message")
      sink(type = "output")
      close(main_log)
    },
    error = function(e) {
      vebcat("Error occurred in node", j, ":", e$message, color = "nonFatalError")
      sink(type = "message")
      sink(type = "output")
      close(main_log)
    }
  )

  rm(list = setdiff(ls(), "j"))

  invisible(gc())

  list(
    iteration = j
  )
}
