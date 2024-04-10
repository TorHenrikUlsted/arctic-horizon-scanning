setup_climate <- function(shapefile, iteration, plot.show, verbose, warn.file, err.file) {
  vebcat("Initiating data acquisition protocol", color = "funInit")
  
  withCallingHandlers(
    {
      biovars_world <- get_wc_data(
        var = climate.var, 
        res = climate.res, 
        plot.show = plot.show, 
        verbose = verbose
      )
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
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when scaling biovars_world", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when scaling biovars_world", iteration = iteration)
  )
  
  vebcat("Acquiring biovars_region.", veb = verbose)
  
  withCallingHandlers(
    {
      biovars_region <- wc_to_region(
        biovars_world,
        shapefile = shapefile, 
        projection = "longlat",
        plot.show = plot.show,
        verbose = verbose
      )
    },
    warning =  function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when acquiring region data", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when acquiring region data", iteration = iteration)
  )
  
  coord_uncertainty <- calc_coord_uncertainty(
    region = biovars_region,
    unit.out = "m",
    dir.out = "./outputs/hypervolume/data_acquisition/logs",
    verbose = verbose
  )
  
  vebcat("Data acquisition protocol completed successfully", color = "funSuccess")
  
  return(list(
    world = biovars_world,
    region = biovars_region
  ))
}

setup_hv_region <- function(biovars_region, out.dir, method) {
  vebcat("Setting up region Hypervolume", color = "funInit")

  region_out <- paste0(paste0(out.dir), "/hypervolume-", method, ".rds")

  if (file.exists(region_out)) {
    catn("Region hypervolume already exists, reading file.")
    region_hv <- readRDS(region_out)
  } else {
    log_dir <- "./outputs/setup/region/logs"
    create_dir_if(log_dir)

    region_log_out <- paste0(log_dir, "/region-hv-", method, "-output.txt")
    region_log_msg <- paste0(log_dir, "/region-hv-", method, "-message.txt")
    create_file_if(c(region_log_out, region_log_msg))

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

setup_hv_sequence <- function(hv.method, hv.accuracy, hv.incl.t, verbose = TRUE) {
  vebcat("Initiating hypervolume sequence setup.", color = "funInit")

  sp_list_setup <- list.files("./outputs/filter/test/test-small/chunk/species", full.names = TRUE)
  
  hv_setup_dir <- "./outputs/setup/hypervolume"
  create_dir_if(hv_setup_dir)

  setup_check <- paste0(hv_setup_dir, "/logs/setup-check.txt")
  low_file <- paste0(hv_setup_dir, "/logs/peak-mem-low.txt")
  high_file <- paste0(hv_setup_dir, "/logs/peak-mem-high.txt")

  min_disk_space <- get_disk_space("/export", units = "GB") * 0.2

  if (file.exists(low_file)) {
    catn("Low peak ram setup already run.")
  } else {
    catn("Running low peak ram setup by using", highcat(sp_list_setup[[1]]), "wait time: 3 min.")

    create_file_if(setup_check)

    ram_control <- start_mem_tracking(file.out = setup_check, stop_file = paste0(hv_setup_dir, "/stop-file.txt"))

    hypercolume_sequence(
      spec.list = sp_list_setup[[1]],
      method = hv.method, # box approx 13 min, gaussian 1 hours 10 minutes
      accuracy = hv.accuracy,
      hv.projection = "longlat",
      proj.incl.t = hv.incl.t,
      iterations = 1,
      min.disk.space = min_disk_space,
      hv.dir = hv_setup_dir,
    )

    stop_mem_tracking(ram_control, low_file, paste0(hv_setup_dir, "/stop-file.txt"))
  }

  if (file.exists(high_file)) {
    catn("high peak ram setup already run.")
  } else {
    catn("Running high peak ram setup using", highcat(sp_list_setup[[2]]), "wait time: 25 min.")
    create_file_if(setup_check)

    ram_control <- start_mem_tracking(file.out = setup_check, stop_file = paste0(hv_setup_dir, "/stop-file.txt"))

    # Run a hypervolume sequence of sax. opp.

    hypercolume_sequence(
      spec.list = sp_list_setup[[3]], # list of strings
      method = hv.method, # box approx 13 min, gaussian 1 hours 10 minutes
      accuracy = hv.accuracy,
      hv.projection = "longlat",
      proj.incl.t = hv.incl.t,
      iterations = 1,
      min.disk.space = min_disk_space,
      hv.dir = hv_setup_dir,
    )

    stop_mem_tracking(ram_control, high_file, paste0(hv_setup_dir, "/stop-file.txt"))
  }

  peak_mem_low <- as.numeric(readLines(low_file))
  peak_mem_high <- as.numeric(readLines(high_file))

  rm(sp_list_setup)
  invisible(gc())

  vebcat("hypervolume sequence setup completed successfully.", color = "funSuccess")

  return(list(
    low = peak_mem_low,
    high = peak_mem_high
  ))
}