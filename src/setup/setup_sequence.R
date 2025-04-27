source_all("./src/setup/components")
if (config$simulation$example) {
  source_all("./example/src/setup")
} else {
  source_all("./src/setup/custom_setup")
}

setup_sequence <- function(approach = "precautionary", coord.uncertainty = NULL, hv.method, hv.accuracy, hv.incl.threshold, hv.dims = NULL, cores.max = 1, force.seq = FALSE, verbose = FALSE) {
  setup_dir <- "./outputs/setup"
  setup_wrangle_dir <- paste0(setup_dir, "/wrangle")
  setup_log <- paste0(setup_dir, "/logs")
  seq_set_file <- paste0(setup_log, "/setup-completed.txt")
  
  if (!is.null(force.seq) && ("all" %in% force.seq || "setup" %in% force.seq)) {
    catn("Forcing Setup sequence.")
    if (file.exists(seq_set_file)) file.remove(seq_set_file)
  }

  if (!file.exists(seq_set_file)) {
    vebcat("Initiating Setup Sequence", color = "seqInit")
    
    mdwrite(
      config$files$post_seq_md,
      text = "1;Setup Sequence"
    )
    
  }

    wfo_speed <- check_system_speed(
      df.path = "./resources/data-raw/test/speed-test-species.csv",
      test.name = "wfo-speed",
      sample.size = NULL,
      cores.max = 1,
      fun = wfo_speed,
      verbose = verbose
    )

    system.speed.wfo <<- wfo_speed

    invisible(gc())

    setup_raw_data(
      column = "interimName",
      cores.max = cores.max,
      verbose = verbose,
      counter = 500
    )
    
    
    if (file.exists(seq_set_file)) {
      vebcat("Setup sequence already run.", color = "proSuccess")
    } else {
    warn_out <- paste0(setup_log, "/warning.txt")
    err_out <- paste0(setup_log, "/error.txt")
    create_file_if(warn_out, err_out)
    
    save_dir <- build_climate_path()
    bw_out <- paste0(save_dir, "/biovars-world-subset.tif")
    br_out <- paste0(save_dir, "/biovars-region-subset.tif")

    sp_dir <- filter_sequence(
      spec.known = list(
        name = "test_known"
      ),
      spec.unknown = list(
        name = "test_small"
      ),
      approach = approach,
      coord.uncertainty = coord.uncertainty,
      cores.max = cores.max,
      force.seq = force.seq,
      verbose = verbose
    )

    rm(sp_dir)
    invisible(gc())

    region_shape <- setup_region()

    biovars <- setup_climate(
      region_shape,
      iteration = 1,
      show.plot = FALSE,
      verbose = verbose,
      warn.file = warn_out,
      err.file = err_out
    )

    if (is.null(hv.dims)) {
      analyzed_data <- analyze_correlation(
        biovars$region,
        dir.out = paste0(setup_dir, "/correlation"),
        verbose = verbose
      )

      vebcat("Check the correlation matrix and pick climate variables, Stopping process", color = "fatalError")
      stop("Choose dimensions based on the correlation matrix")
    } else {
      catn("loading biovars.")

      biovars$world <- terra::subset(biovars$world, hv.dims)
      if (!file.exists(bw_out)) writeRaster(biovars$world, bw_out)

      biovars$region <- terra::subset(biovars$region, hv.dims)
      if (!file.exists(br_out)) writeRaster(biovars$region, br_out)

      region_hv <- setup_hv_region(
        biovars$region,
        out.dir = build_climate_path(),
        method = hv.method
      )

      rm(biovars)
      invisible(gc())

      setup_hv_sequence(
        hv.method = hv.method,
        hv.accuracy = hv.accuracy,
        hv.dims = hv.dims,
        hv.incl.threshold = hv.incl.threshold,
        verbose = verbose
      )
    }

    create_file_if(seq_set_file)

    vebcat("Setup Sequence Completed Successfully", color = "seqSuccess")
  }
}
