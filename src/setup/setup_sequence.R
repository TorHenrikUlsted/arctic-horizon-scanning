source_all("./src/setup/components")
source_all("./src/setup/custom_setup")

setup_sequence <- function(hv.method, hv.accuracy, hv.incl.t, hv.dims = NULL, cores.max = 1, verbose = FALSE) {
  
  vebcat("Initiating Setup Sequence", color = "seqInit")
  
  setup_dir <- "./outputs/setup"
  setup_log <- paste0(setup_dir, "/log")
  
  warn_out <- paste0(setup_log, "/warning.txt")
  err_out <- paste0(setup_log, "/error.txt")
  
  create_file_if(c(warn_out, err_out))
  
  wfo_speed <- check_system_speed(
    df.path = "./resources/data-raw/speed-test-species.csv",
    test.name = "wfo-speed",
    sample.size = NULL,
    cores.max = 1,
    fun = wfo_speed,
    verbose = verbose
  )
  
  system.speed.wfo <<- wfo_speed
  
  setup_raw_data(
    column = "rawName",
    cores.max = cores.max,
    verbose = verbose,
    counter = 10
  )
  
  region_shape <- setup_region()
  
  biovars <- data_acquisition(
    region_shape,
    iteration = 1, 
    show.plot = FALSE, 
    verbose = verbose, 
    warn.file = warn_out,
    err.file = err_out
  )
  
  analyzed_data <- analyze_correlation(
    biovars$region, 
    file.out = paste0(setup_dir, "/correlation"), 
    verbose = verbose
  )
  
  if (is.null(hv.dims)) {
    vebcat("Check the correlation matrix and pick climate variables, Stopping process", color = "fatalError")
    stop("")
  } else {
    bw_out <- paste0(setup_dir, "/region/biovars-world-subset.tif")
    br_out <- paste0(setup_dir, "/region/biovars-region-subset.tif")
    
    biovars_world <- terra::subset(biovars$world, hv.dims)
    if (!file.exists(bw_out)) writeRaster(biovars_world, bw_out)
    
    biovars_region <- terra::subset(biovars$region, hv.dims)
    if (!file.exists(br_out)) writeRaster(biovars_region, br_out)
    
    region_hv <- setup_region_hv(
      biovars$region, 
      out.dir = paste0(setup_dir, "/region"), 
      method = hv.method
    )
    
    setup_hv_sequence(
      hv.method = hv.method, 
      hv.accuracy = hv.accuracy, 
      hv.incl.t = hv.incl.t,
      verbose = verbose
    )
  }
  
  vebcat("Setup Sequence Completed Successfully", color = "seqSuccess")
}
