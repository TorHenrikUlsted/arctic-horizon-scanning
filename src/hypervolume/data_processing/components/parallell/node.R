node_processing <- function(j, sp_list, proj.incl.t, method, projection, show.plot, verbose, min.disk.space) {
  current_disk_space <- get_disk_space("/home", units = "GB")

  warn <- function(w, warn_txt) {
    warn_msg <- conditionMessage(w)
    sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
    cat(warn_txt, j, ":", warn_msg, "\n")
    sink()
    invokeRestart(findRestart("muffleWarning"))
  }

  err <- function(e, err_txt) {
    err_msg <- conditionMessage(e)
    sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
    cat(err_txt, j, ":", err_msg, "\n")
    sink()
    stop(e)
  }

  # Check if the current disk space is less than the minimum required
  if (current_disk_space <= min.disk.space) stop("Insufficient disk space. Stopping processing.")

  spec.name <- names(sp_list)[j]

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

  acq_data <- data_acquisition(show.plot, method, verbose = verbose, iteration = j, warn, err)

  invisible(gc())

  ana_data <- data_analysis(biovars_world = acq_data[[1]], biovars_region = acq_data[[2]], biovars_floreg = acq_data[[3]], method, show.plot, verbose = verbose, iteration = j, warn, err)

  invisible(gc())

  processed_data <- data_processing(sp_list[[j]], biovars_world = ana_data[[1]], spec.name, method, projection, verbose = verbose, iteration = j, warn, err)
  invisible(gc())

  analyzed_hv <- hv_analysis(processed_data, biovars_region = ana_data[[2]], region_hv = ana_data[[4]], method, spec.name, proj.incl.t, verbose = verbose, iteration = j, warn, err)

  final_res <- data.frame(
    Species = spec.name,
    observations = analyzed_hv[[1]],
    dimensions = analyzed_hv[[2]],
    samplesPerPoint = analyzed_hv[[3]],
    randomPoints = analyzed_hv[[4]],
    excluded = analyzed_hv[[5]],
    jaccard = analyzed_hv[[6]][[1]],
    sørensen = analyzed_hv[[6]][[2]],
    volumeSpecies = analyzed_hv[[6]][[3]],
    volumeRegion = analyzed_hv[[6]][[4]],
    overlapRegionBox = 1 - analyzed_hv[[6]][[4]],
    includedOverlapT1 = sum(analyzed_hv[[7]] == T) / length(analyzed_hv[[7]]),
    includedOverlapT25 = sum(analyzed_hv[[8]] == T) / length(analyzed_hv[[8]]),
    includedOverlapT5 = sum(analyzed_hv[[9]] == T) / length(analyzed_hv[[9]])
  )
  # If we want the number of true and false values we can do: round(observations * includedOverlap) = true values, and: observations - true values = false values

  fwrite(final_res, paste0("./outputs/hypervolume/stats/", method, "-stats.csv"), append = T, bom = T)

  end_timer(node_timer)

  sink(type = "message")
  sink(type = "output")
  close(log_file)

  rm(list = setdiff(ls(), "j"))

  invisible(gc())

  list(
    iteration = j,
    warning = FALSE,
    error = FALSE
  )
}
