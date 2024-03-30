setup_region_hv <- function(biovars_region, out.dir, name, method) {
  vebcat("Setting up region Hypervolume", color = "funInit")
  
  region_filename <- paste0("./outputs/setup/region/", name, "/hypervolume-", method,".rds")
  
  if (file.exists(region_filename)) {
    catn("Region hypervolume already exists, reading file.")
    region_hv <- readRDS(region_filename)
  } else {
    create_dir_if("./outputs/setup/region/logs")
    create_dir_if(paste0("./outputs/setup/region/", name))
    region_log_out <- paste0("./outputs/setup/region/logs/", name, "-", method, "-output.txt")
    
    create_file_if(region_log_out)
    
    region_log_msg <- paste0("./outputs/setup/region/logs/", name, "-", method, "-message.txt")
    create_file_if(region_log_msg)
    
    try(region_log_out <- file(region_log_out, open = "at"))
    try(region_log_msg <- file(region_log_msg, open = "at"))
    sink(region_log_out, type = "output")
    sink(region_log_msg, type = "message")
    
    region_hv_timer <- start_timer("region_hv_timer")
    
    region_hv <- analyze_region_hv(biovars_region, out.dir, name, method = method, verbose = T)
    
    end_timer(region_hv_timer)
    
    sink(type = "message")
    sink(type = "output")
    close(region_log_out)
    close(region_log_msg)
  }
  
  return(region_hv)
}

setup_hv_sequence <- function(min_disk_space, verbose = T) {
  vebcat("Initiating hypervolume sequence setup.", color = "funInit")
  
  sp_list_setup <- list.files("./outputs/filter/test/test-small/chunk/species", full.names = TRUE)
  hv_setup_dir <- "./outputs/setup/hypervolume"
  create_dir_if(hv_setup_dir)
  setup_check <- paste0(hv_setup_dir, "/logs/setup-check.txt")
  low_file <- paste0(hv_setup_dir, "/logs/peak-mem-low.txt")
  high_file <- paste0(hv_setup_dir, "/logs/peak-mem-high.txt")
  
  
  
  if (file.exists(low_file)) {
    catn("Low peak ram setup already run.")
  } else {
    catn("Running low peak ram setup by using", highcat(sp_list_setup[[1]]), "wait time: 3 min.")
    
    create_file_if(setup_check)
    
    ram_control <- start_mem_tracking(file.out = setup_check, stop_file = paste0(hv_setup_dir, "/stop-file.txt"))
    
    parallell_processing(
      spec.list = sp_list_setup[[1]],
      method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
      accuracy = "accurate",
      hv.projection = "laea",
      proj.incl.t = 0.5,
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
    
    #Run a hypervolume sequence of sax. opp.
    
    parallell_processing(
      spec.list = sp_list_setup[[3]], # list of strings
      method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
      accuracy = "accurate",
      hv.projection = "laea",
      proj.incl.t = 0.5,
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