parallel_spec_dirs <- function(spec.dirs, dir, hv.project.method, fun, verbose = FALSE) {
  cat("Setting up parallel function.\n")
  
  ram_log <- paste0(dir, "/ram-usage.txt")
  create_file_if(ram_log)
  
  node_dir <- paste0(dir, "/nodes")
  create_dir_if(node_dir)
  
  # Get max core usage
  used_mem <- get_mem_usage(type = "used", format = "gb")
  remain_mem <- ((mem_limit / 1024^3) - used_mem)
  
  check_mem_need <- function()  {
    mem_need <- 5
    
    return(mem_need)
  }
  
  mem_per_core_gb <- check_mem_need()
  
  cores_max <- floor(remain_mem / mem_per_core_gb)
  
  cores_max <- min(length(spec.dirs), cores_max, total_cores)
  
  cat("Creating cluster of cores. \n")
  
  cl <- makeCluster(cores_max)
  
  clusterExport(cl, c("spec.dirs", "fun", "hv.project.method", "node_dir", "ram_log"), envir = environment())
  
  clusterEvalQ(cl, {
    source("./src/utils/utils.R")
    source("./src/visualize/visualize.R")
  })
  
  current_disk_space <- get_disk_space("/export", units = "GB")
  
  # If iterations is not provided, start from the highest saved iteration
  i <- 1
  end <- length(spec.dirs)
  iterations <- i:end
  
  ram_msg <- FALSE
  # RAM check
  mem_used_gb <- get_mem_usage(type = "used", format = "gb")
  mem_limit_gb <- mem_limit / 1024^3
  
  cat("Running parallel with", cores_max, "cores.\n")
  cat("Max iterations is", end,"\n")
  
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
    
    node_con <- file(node_log, open = "at")
    if (!inherits(node_con, "connection")) {
      stop("Failed to open log file: ", node_log)
    }
    
    sink(node_con, type = "output")
    sink(node_con, type = "message")
    
    sp_dir <- spec.dirs[[i]]
    sp_name <- basename(sp_dir)
    sp_filename <- paste0(sp_dir, "/", hv.project.method, ".tif")
    
    cat("Running node for:\n")
    cat("Iteration",i,"\n")
    cat(sp_name, "\n")
    cat(sp_dir, "\n")
    cat(sp_filename, "\n")
    
    dt <- fun(spec.filename = sp_filename)
    
    sink(type = "message")
    sink(type = "output")
    close(node_con)
    
    return(data = dt)
  })
  
  stopCluster(cl)
  
  return(res)
}

parallel_spec_handler <- function(spec.dirs, hv.project.method = "inclusion-0.5", sub.dir, n = NULL, fun, verbose = FALSE) {
  
  cat(blue("Acquiring", basename(sub.dir), "values for all", hv.project.method, "rasters.\n"))
  
  create_dir_if(sub.dir)
  
  values_sub_dir <- paste0(sub.dir, "/", hv.project.method)
  create_dir_if(values_sub_dir)
  
  values_log <- paste0(values_sub_dir, "/", basename(sub.dir), ".csv")
  
  if (verbose) {
    print(sub.dir)
    print(basename(sub.dir))
    print(values_sub_dir)
    print(values_log)
  }
  
  if (file.exists(values_log)) {
    cat("Found coverage value table:", values_log, "\n")
    
    combined_values <- fread(values_log)
    
  } else {
    cat("No previous table found, acquiring coverage values for", hv.project.method, "method.\n")
    
    cat("Setting up function.\n")
    
    combined_values <- data.table()
    
    values <- parallel_spec_dirs(spec.dirs, values_sub_dir, hv.project.method, fun = fun, verbose = verbose)
    
    cat("Finishing up \n")
    
    combined_values <- rbindlist(values)
    
    combined_values <- combined_values[order(-combined_values$coverage), ]
    
    fwrite(combined_values, values_log, bom = TRUE)
  }
  
  if (is.null(n)) {
    cat("n is not used, returning the whole table.\n")
  } else {
    combined_values <- combined_values[1:n, ]
  }
  
  cat(cc$lightGreen("Successfully acquired coverage values for all probability rasters.\n"))
  
  return(combined_values)
}