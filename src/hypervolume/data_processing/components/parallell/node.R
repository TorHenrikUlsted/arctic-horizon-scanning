node_processing <- function(j, sp_list, method, show_plot, verbose) {
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

  acq_data <- data_acquisition(show_plot, verbose, iteration = j)

  invisible(gc())

  ana_data <- data_analysis(biovars_world = acq_data[[1]], biovars_region = acq_data[[2]], biovars_floreg = acq_data[[3]], iteration = j)

  invisible(gc())

  processed_data <- data_processing(sp_list[[j]], biovars_world = ana_data[[1]], projection, verbose, iteration = j)

  invisible(gc())

  analyzed_hv <- hv_analysis(processed_data, biovars_region = ana_data[[2]], region_hv = ana_data[[4]], method, spec.name, verbose, iteration = j)


  final_res <- data.frame(
    Species = spec.name,
    jaccard = analyzed_hv[[2]][[1]],
    sørensen = analyzed_hv[[2]][[2]],
    volumeSpecies = analyzed_hv[[2]][[3]],
    volumeRegion = analyzed_hv[[2]][[4]],
    overlapRegion = 1 - analyzed_hv[[2]][[4]],
    randomPoints = analyzed_hv[[3]],
    samplesPerPoint = analyzed_hv[[4]],
    observations = analyzed_hv[[5]],
    dimensions = analyzed_hv[[6]],
    maxInclusionPoints = length(analyzed_hv[[1]]),
    includedPoints = sum(analyzed_hv[[1]] == T),
    excludedPoints = sum(analyzed_hv[[1]] == F),
    includedOverlap = sum(analyzed_hv[[1]] == T) / length(analyzed_hv[[1]])
  )
  
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
