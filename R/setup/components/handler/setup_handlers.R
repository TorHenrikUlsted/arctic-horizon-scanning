setup_climate <- function(shapefile, iteration, show.plot = FALSE, verbose = FALSE, warn.file, err.file) {
  vebcat("Initiating climate setup protocol", color = "funInit")

  withCallingHandlers(
    {
      biovars_world <- load_climate_data(
        database = config$simulation$climate,
        show.plot = show.plot,
        verbose = verbose
      )
      invisible(gc())
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when getting worldClim data", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when getting worldClim data", iteration = iteration)
  )

  vebcat("Scaling biovars", veb = verbose)

  withCallingHandlers(
    {
      biovars_world <- scale_biovars(
        biovars_world,
        verbose = verbose
      )
      invisible(gc())
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when scaling biovars_world", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when scaling biovars_world", iteration = iteration)
  )

  vebcat("Acquiring biovars_region.", veb = verbose)

  withCallingHandlers(
    {
      biovars_region <- climate_to_region(
        biovars_world,
        shapefile = shapefile,
        projection = "longlat",
        show.plot = show.plot,
        verbose = verbose
      )
      invisible(gc())
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when acquiring region data", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when acquiring region data", iteration = iteration)
  )

  withCallingHandlers(
    {
      if (!is.null(config$projection$raster_scale_m)) {
        coord_uncertainty <- config$projection$raster_scale_m
      } else {
        coord_uncertainty <- calc_coord_uncertainty(
          region = biovars_region,
          projection = config$simulation$projection,
          unit.out = "m",
          dir.out = build_climate_path(),
          verbose = verbose
        )
      }
      
      invisible(gc())
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when calculating coordinate uncertainty", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when calculating coordinate uncertainty", iteration = iteration)
  )

  vebcat("Climate setup protocol completed successfully", color = "funSuccess")

  return(list(
    world = biovars_world,
    region = biovars_region
  ))
}

setup_hv_region <- function(biovars_region, out.dir, method) {
  on.exit({
    rm(list = ls(environment()))
    gc(full = TRUE)
  })
  
  vebcat("Setting up region Hypervolume", color = "funInit")

  region_out <- paste0(paste0(out.dir), "/hypervolume-", method, ".rds")

  if (file.exists(region_out)) {
    catn("Region hypervolume already exists, reading file.")
    region_hv <- readRDS(region_out)
  } else {
    log_dir <- paste0(dirname(out.dir), "/logs")
    create_dir_if(log_dir)

    region_log_out <- paste0(log_dir, "/region-hv-", method, "-output.txt")
    region_log_msg <- paste0(log_dir, "/region-hv-", method, "-message.txt")
    create_file_if(region_log_out, region_log_msg)

    try(region_log_out <- file(region_log_out, open = "at"))
    try(region_log_msg <- file(region_log_msg, open = "at"))
    sink(region_log_out, type = "output")
    sink(region_log_msg, type = "message")

    region_hv_timer <- start_timer("region_hv_timer")

    region_hv <- analyze_region_hv(
      biovars_region,
      file.out = region_out,
      method = method,
      verbose = TRUE
    )

    end_timer(region_hv_timer)

    sink(type = "message")
    sink(type = "output")
    close(region_log_out)
    close(region_log_msg)
  }

  return(region_hv)
}

setup_hv_sequence <- function(hv.method, hv.accuracy, hv.dims, hv.incl.threshold = 0.5, verbose = TRUE) {
  vebcat("Initiating hypervolume sequence setup.", color = "funInit")

  sp_list_setup <- list.files("./outputs/filter/test-small/chunk/species", full.names = TRUE)

  hv_setup_dir <- "./outputs/setup/hypervolume"
  hv_logs_dir <- paste0(hv_setup_dir, "/", hv.method, "-sequence/setup-check")
  create_dir_if(hv_logs_dir) # Recursive

  stop_file <- paste0(hv_logs_dir, "/stop-file.txt")
  low_file <- paste0(hv_logs_dir, "/peak-mem-low.txt")
  high_file <- paste0(hv_logs_dir, "/peak-mem-high.txt")
  time_low_file <- paste0(hv_logs_dir, "/time-low.txt")
  time_high_file <- paste0(hv_logs_dir, "/time-high.txt")
  avg_time_file <- paste0(hv_logs_dir, "/avg-time.txt")

  min_disk_space <- get_disk_space("/export", units = "GB") * 0.2

  if (file.exists(low_file)) {
    catn("Low peak ram setup already run.")
  } else {
    catn("Running low peak ram setup by using", highcat(sp_list_setup[[1]]), "wait time: 3 min.")

    tryCatch({
      catn("Starting memory tracker")
      ram_control <- start_mem_tracking(low_file, stop_file)
      time_low <- start_timer("low_mem_test")

      hypervolume_sequence(
        spec.list = sp_list_setup,
        iterations = "Abelmoschus manihot",
        min.disk.space = min_disk_space,
        verbose = TRUE,
        hv.dir = hv_setup_dir,
        hv.method = hv.method,
        hv.accuracy = hv.accuracy,
        hv.dims = hv.dims,
        hv.incl.threshold = hv.incl.threshold
      )
    }, error = function(e) {
      vebcat("Error when running low memory hypervolume sequence", color = "fatalError")
      stop(e)
    }, finally = {
      stop_mem_tracking(ram_control, stop_file)
      time_low_res <- end_timer(time_low)
      writeLines(as.character(time_low_res), time_low_file)
    })
  }

  if (file.exists(high_file)) {
    catn("high peak ram setup already run.")
  } else {
    catn("Running high peak ram setup using", highcat(sp_list_setup[[3]]), "wait time: 25 min.")
    create_file_if(high_file)

    tryCatch({
      ram_control <- start_mem_tracking(high_file, stop_file)
      time_high <- start_timer("high_mem_test")

      # Run a hypervolume sequence of sax. opp.
      hypervolume_sequence(
        spec.list = sp_list_setup,
        iterations = "Saxifraga oppositifolia",
        min.disk.space = min_disk_space,
        verbose = TRUE,
        hv.dir = hv_setup_dir,
        hv.method = hv.method,
        hv.accuracy = hv.accuracy,
        hv.dims = hv.dims,
        hv.incl.threshold = hv.incl.threshold
      )
    }, error = function(e) {
      vebcat("Error when running high memory hypervolume sequence", color = "fatalError")
      stop(e)
    }, finally = {
      stop_mem_tracking(ram_control, stop_file)
      time_high_res <- end_timer(time_high)
      writeLines(as.character(time_high_res), time_high_file)
    })
  }

  peak_mem_low <- as.numeric(readLines(low_file))
  peak_mem_high <- as.numeric(readLines(high_file))
  
  if(file.exists(time_low_file) && file.exists(time_high_file)) {
    time_low <- as.numeric(readLines(time_low_file))
    time_high <- as.numeric(readLines(time_high_file))
    avg_time <- (time_low + time_high) / 2
    
    # Convert to minutes for more readable display
    time_low_min <- time_low / 60
    time_high_min <- time_high / 60 
    avg_time_min <- avg_time / 60
    
    writeLines(as.character(round(avg_time, 2)), avg_time_file)
    
    catn("\nTiming Results:")
    catn("Low Memory Species Time:", highcat(round(time_low_min, 2)), "minutes")
    catn("High Memory Species Time:", highcat(round(time_high_min, 2)), "minutes") 
    catn("Average Processing Time:", highcat(round(avg_time_min, 2)), "minutes")
  }

  rm(sp_list_setup)
  invisible(gc())

  vebcat("hypervolume sequence setup completed successfully.", color = "funSuccess")

  return(list(
    low = peak_mem_low,
    high = peak_mem_high,
    time = list(
      low = if(file.exists(time_low_file)) as.numeric(readLines(time_low_file)) else NULL,
      high = if(file.exists(time_high_file)) as.numeric(readLines(time_high_file)) else NULL,
      avg = if(file.exists(avg_time_file)) as.numeric(readLines(avg_time_file)) else NULL
    )
  ))
}
