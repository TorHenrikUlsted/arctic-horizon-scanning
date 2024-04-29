parallel_spec_dirs <- function(spec.dirs, dir, shape, extra, hv.project.method, fun.init, fun.execute, batch = FALSE, test = 0, node.log.append = TRUE, verbose = FALSE) {

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
    shape = shape,
    extra = extra,
    fun.init = fun.init,
    file.log = peak_log, 
    file.stop = stop_file,
    verbose = verbose
  )
  
  if (length(sp_dirs) == 0) {
    catn(highcat(length(sp_dirs)), "species left after the initiation, returning.")
    return(list(init.res = init_res, exec.res = NULL))
  }
  
  timer_res <- end_timer(init_timer)

  peak_ram <- as.numeric(readLines(peak_log))
  
  max_cores <- calc_num_cores(
    ram.high = peak_ram, 
    verbose =  FALSE
  )
  max_cores <- max_cores$high
  
  if (test > 0) { # use test numbers if above 0
    max_cores <-  min(test, max_cores, total_cores)
  } else {
    max_cores <- min(length(sp_dirs), max_cores, total_cores)
  }
  
  catn("Creating cluster of", max_cores, "core(s).")
  
  cl <- makeCluster(max_cores)
  
  clusterEvalQ(cl, {
    source("./src/utils/utils.R")
    source("./src/setup/setup_sequence.R")
    source_all("./src/visualize/components")
  })
  
  # Get the exports that are not null
  export_vars <- c("sp_dirs", "shape", "hv.project.method", "fun.execute", "batch", "node.log.append", "test", "verbose" , "node_dir", "ram_log", "extra")
  
  #export_vars <- export_vars[sapply(export_vars, function(x) !is.null(get(x)))]
  
  vebcat("export_vars:", veb = verbose)
  vebprint(export_vars, veb = verbose)
  
  clusterExport(cl, export_vars, envir = environment())

  current_disk_space <- get_disk_space("/export", units = "GB")
  
  # Check for batched approach
  if (batch) {
    catn("Splitting species directories into batches.")
    # Calculate the number of species per batch
    sp_per_core <- ceiling(length(sp_dirs) / max_cores)
    
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
  calculate_etc(timer_res, max_cores, length(spec.dirs))
  
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
      
      vebprint(shape, verbose, "Shape:")
      if (!is.null(shape)) {
        region <- load_region(shape)
        region <- handle_region(region)
      } else {
        region <- NULL
      }
      
      vebprint(extra, verbose, "Extra:")
      
      catn("\nParallel processing execute function.")
      dt <- fun.execute(
        spec.filename = sp_filename, 
        region = region,
        extra = extra,
        verbose = verbose
      )
      
      catn("Parallel execute function completed successfully.")
      
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
    return(list(
      init.res = init_res,
      exec.res = res
    ))
  }
}

parallel_spec_handler <- function(spec.dirs, dir, shape = NULL, extra = NULL, hv.project.method = "0.5-inclusion", n = NULL, out.order = NULL, fun, batch = FALSE, test = 0, node.log.append = TRUE, verbose = FALSE) {
  if (verbose) {
    print_function_args()
  }
  input_functions <- fun()
  
  if (!exists("execute", where = input_functions)) {
    vebcat("Missing execute function to run parallel process.", color = "fatalError")
    stop()
  } else {
    vebcat("Found execute function.")
    fun_execute <- input_functions$execute
  }
  
  if (!exists("init", where = input_functions)) {
    vebcat("Setting init function be the same as execute.")
    fun_init <- input_functions$execute
  } else {
    vebcat("Found init function.")
    fun_init <- input_functions$init
  } 
    
  if (!exists("process", where = input_functions)) {
    vebcat("Setting process function to use single function.")
    fun_process <- parallel_process_single
  } else {
    vebcat("Found process function.")
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
      shape = shape,
      extra = extra,
      hv.project.method = hv.project.method, 
      fun.init = fun_init,
      fun.execute = fun_execute,
      batch = batch,
      test = test,
      node.log.append = node.log.append,
      verbose = verbose
    )
    
    catn("Length of init results:", highcat(length(parallel_res$init.res)))
    catn("Length of execute results:", highcat(length(parallel_res$exec.res)))
    
    vebprint(parallel_res$init.res, verbose, "Parallel init result:")
    vebprint(head(parallel_res$exec.res, 1), verbose, "Parallel execute head(result, 1):")
    
    catn("Running post processing.")
    # Run the post process
    if (is.null(parallel_res$exec.res)) {
      process_res <- parallel_res$init.res
    } else {
      process_res <- fun_process(
        parallel.res = parallel_res,
        extra = extra,
        verbose = verbose
      )  
    }
    
    vebprint(head(process_res, 3), verbose, "Processed data:")
    
    fwrite(process_res, values_log, bom = TRUE)
  }
  
  if (!is.null(out.order)) {
    process_res <- process_res[order(-process_res[[out.order]]), ]
  }
  
  if (is.null(n)) {
    catn("n is not used, returning the whole table")
  } else {
    process_res <- process_res[1:n, ]
  }
  
  vebcat("Successfully acquired", basename(dir), "for all",  hv.project.method, "rasters.", color = "funSuccess")
  
  return(process_res)
}

parallel_init <- function(spec.filename, shape, extra, hv.project.method, fun.init, file.log, file.stop, verbose) {
 
  tryCatch({
    # Initiate memory control
    ram_control <- start_mem_tracking(file.out = file.log, file.stop = file.stop)
    
    if (!is.null(shape)) {
      region <- load_region(shape)
      region <- handle_region(region)
    } else {
      region <- NULL
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
  catn("Using the single processor.")
  merged_dt <- rbindlist(list(parallel.res$init.res, rbindlist(parallel.res$exec.res)), fill = TRUE)
  
  return(merged_dt)
}