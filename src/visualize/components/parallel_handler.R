parallel_spec_dirs <- function(spec.dirs, dir, region, hv.project.method, fun.init, fun.execute, batch = FALSE, test = 0, node.log.append = TRUE, verbose = FALSE) {

  if (verbose) {
    print_function_args()
  }
  
  catn("Setting up parallel function")
  
  ram_log <- paste0(dir, "/ram-usage.txt")
  peak_log <- paste0(dir, "/ram_peak.txt")
  stop_file <- paste0(dir, "/stop-file.txt")
  
  create_file_if(c(ram_log, peak_log))
  
  node_dir <- paste0(dir, "/nodes")
  create_dir_if(node_dir)
  
  # Get max core usage
  used_mem <- get_mem_usage(type = "used", format = "gb")
  remain_mem <- ((mem_limit / 1024^3) - used_mem)
  
  # Set up names for the initiation
  sp_dir <- spec.dirs[1]
  sp_dirs <- spec.dirs[-1]
  sp_filename <- paste0(sp_dir, "/", hv.project.method, ".tif")
  
  catn("Processing init function.")
  init_timer <- start_timer("init timer")
  
  init_res <- parallel_init(
    spec.filename = sp_filename, 
    region = region, 
    fun.init = fun.init,
    file.log = peak_log, 
    file.stop = stop_file,
    verbose = verbose
  )
  
  timer_res <- end_timer(init_timer)

  peak_mem_core <- as.numeric(readLines(peak_log))
  
  cores_max <- floor(remain_mem / peak_mem_core)
  
  if (test > 0) { # use test numbers if above 0
    cores_max <-  min(test, cores_max, total_cores)
  } else {
    cores_max <- min(length(sp_dirs), cores_max, total_cores)
  }
  
  catn("Creating cluster of", cores_max, "core(s).")
  
  cl <- makeCluster(cores_max)
  
  clusterEvalQ(cl, {
    source("./src/utils/utils.R")
    source_all("./src/visualize/components")
    source("./src/setup/components/region/import_regions.R")
  })
  
  # Get the exports that are not null
  export_vars <- c("sp_dirs", "region", "hv.project.method", "fun.execute", "batch", "node.log.append", "test", "verbose" , "node_dir", "ram_log")
  
  export_vars <- export_vars[sapply(export_vars, function(x) !is.null(get(x)))]
  
  vebcat("export_vars:", veb = verbose)
  vebprint(export_vars, veb = verbose)
  
  clusterExport(cl, export_vars, envir = environment())

  current_disk_space <- get_disk_space("/export", units = "GB") #WIP
  
  # Check for batched approach
  if (batch) {
    catn("Splitting species directories into batches.")
    # Calculate the number of species per batch
    sp_per_core <- ceiling(length(sp_dirs) / cores_max)
    
    # Split tasks into batches
    sp_dirs <- split(sp_dirs, ceiling(seq_along(sp_dirs) / sp_per_core))
  }
  
  # If iterations is not provided, start from the highest saved iteration
  i <- 1
  if (test > 0) end <- test else end <- length(sp_dirs)
  iterations <- i:end
  
  ram_msg <- FALSE
  # RAM check
  mem_used_gb <- get_mem_usage(type = "used", format = "gb")
  mem_limit_gb <- mem_limit / 1024^3
  
  catn("Running parallel with", end, "iteration(s).")
  calculate_etc(timer_res, cores_max, length(spec.dirs))
  
  tryCatch({
    res <- clusterApplyLB(cl, iterations, function(i) {
      while (mem_used_gb >= mem_limit_gb) {
        if (!ram_msg) {
          ram_con <- file(ram_log, open = "a")
          writeLines(paste0("RAM usage ", mem_used_gb, " is above the maximum ", mem_limit_gb, " Waiting with node", i), ram_con)
          close(ram_con)
          ram_msg = TRUE
        }
        Sys.sleep(60)  # Wait for 5 seconds before checking again
        Sys.sleep(runif(1, 0, 1)) # Add random seconds between 0 and 1 to apply difference if multiple nodes are in queue
        mem_used_gb <- get_mem_usage(type = "used", format = "gb")
      }
      
      # Setup node log
      node_log <- paste0(node_dir, "/node-", i, ".txt")
      create_file_if(node_log)
      
      if (node.log.append) {
        node_con <- file(node_log, open = "at")
      } else {
        node_con <- file(node_log, open = "wt")
      }
      if (!inherits(node_con, "connection")) {
        stop("Failed to open log file: ", node_log)
      }
      
      sink(node_con, type = "output")
      sink(node_con, type = "message")
      
      sp_dir <- sp_dirs[[i]]
      sp_name <- basename(sp_dir)
      sp_filename <- paste0(sp_dir, "/", hv.project.method, ".tif")
      
      catn("Running node for:")
      
      if (!batch) {
        catn("Iteration", i)
        catn(sp_name)
        catn(sp_dir)
        catn(sp_filename, "\n")
      } else {
        if (test > 0) sp_filename <- sp_filename[1:3]
        catn("Batch", i)
        catn("Number of species:", length(sp_dir))
        catn("Number of filenames:", length(sp_filename), "\n")
      }
      
      catn("\nParallel processing execute function.")
      dt <- fun.execute(
        spec.filename = sp_filename, 
        region = region,
        verbose = verbose
      )
      
      sink(type = "message")
      sink(type = "output")
      close(node_con)
      
      return(data = dt)
    })
  }, error = function(e) {
    vebcat("An error occurred in the parallel process ~ stopping cluster and closing connections.", color = "fatalError")
    stopCluster(cl)
    closeAllConnections()
    stop(e$message)
  })
  
  stopCluster(cl)
  
  if (batch) {
    return(list(
      init.res = init_res,
      exec.res = res
    ))
  } else {
    return(res)
  }
}

parallel_spec_handler <- function(spec.dirs, dir, region = NULL, hv.project.method = "inclusion-0.5", n = NULL, out.order = NULL, fun, batch = FALSE, test = 0, node.log.append = TRUE, verbose = FALSE) {
  if (verbose) {
    print_function_args()
  }
  input_functions <- fun()
  
  if (!exists("execute", where = input_functions)) {
    vebcat("Missing execute function to run parallel process.", color = "fatalError")
    stop()
  } else {
    fun_execute <- input_functions$execute
  }
  
  if (!exists("init", where = input_functions)) {
    fun_init <- input_functions$execute
  } else {
    fun_init <- input_functions$init
  } 
    
  if (!exists("process", where = input_functions)) {
    fun_process <- parallel_process_single
  } else {
    fun_process <- input_functions$process
  }
  
  vebcat("Acquiring", basename(dir), "values for all", hv.project.method, "rasters", color = "funInit")
  
  create_dir_if(dir)
  
  values_sub_dir <- paste0(dir, "/", hv.project.method)
  create_dir_if(values_sub_dir)
  
  values_log <- paste0(values_sub_dir, "/", basename(dir), ".csv")
  
  if (file.exists(values_log)) {
    catn("Found value table:", values_log)
    
    process_res <- fread(values_log)
    
  } else {
    catn("No previous table found, acquiring values for", hv.project.method, "method")
    
    combined_values <- data.table()
    
    parallel_res <- parallel_spec_dirs(
      spec.dirs = spec.dirs,
      dir = values_sub_dir,
      region = region,
      hv.project.method = hv.project.method, 
      fun.init = fun_init,
      fun.execute = fun_execute,
      batch = batch,
      test = test,
      node.log.append = node.log.append,
      verbose = verbose
    )
    
    vebprint(head(parallel_res, 3), verbose)
    
    catn("Running post processing.")
    
    # Run the post process
    process_res <- fun_process(
      parallel.res = parallel_res,
      verbose = verbose
    )
    
    if (!is.null(out.order)) {
      process_res <- process_res[order(-process_res[[out.order]]), ]
    }
    
    fwrite(process_res, values_log, bom = TRUE)
  }
  
  if (is.null(n)) {
    catn("n is not used, returning the whole table")
  } else {
    process_res <- process_res[1:n, ]
  }
  
  vebcat("Successfully acquired", basename(dir), "for all",  hv.project.method, "rasters.", color = "funSuccess")
  
  return(process_res)
}

parallel_init <- function(spec.filename, region, hv.project.method, fun.init, file.log, file.stop, verbose) {
 
  tryCatch({
    # Initiate memory control
    ram_control <- start_mem_tracking(file.out = file.log, file.stop = file.stop)
    
    if (!is.null(region)) {
      regions <- import_regions(region, "./outputs/visualize/logs/region")
      region <- handle_region(regions[[1]])
    }
    
    vebprint(class(fun.init), veb = verbose)
    vebprint(spec.filename, veb = verbose)
    vebprint(region, veb = verbose)
    vebprint(verbose, veb = verbose)
    
    init_res <- fun.init(
      spec.filename = spec.filename, 
      region = region, 
      verbose = verbose
    )
    
    # Stop the memory tracker
    stop_mem_tracking(ram_control, file.stop = file.stop)
  }, error = function(e) {
    vebcat("Error occurred:\n", e$message, color = "nonFatalError")
    vebcat("Stopping RAM check process.", color = "nonFatalError")
    stop_mem_tracking(ram_control, file.stop = file.stop)
  })
  
  return(init_res)
}

parallel_process_single <- function(parallel.res, verbose) {
  merged_dt <- rbindlist(parallel.res, fill = TRUE)
  
  return(merged_dt)
}