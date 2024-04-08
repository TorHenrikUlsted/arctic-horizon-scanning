setup_parallel <- function(par.dir, spec.list, iterations, cores.max, cores.max.high, min.disk.space, verbose = FALSE, custom.exports, custom.evals) {
  parallell_timer <- start_timer("parallell_timer")

  vebcat("Setting up hypervolume files and folders.", veb = verbose)

  logs_dir <- paste0(par.dir, "/logs")

  create_dir_if(c(par.dir, logs_dir, stats_dir))

  ram_usage <- paste0(logs_dir, "/ram-usage.txt")
  node_it_file <- paste0(logs_dir, "/node-iterations.txt")
  highest_it_file <- paste0(logs_dir, "/highest_iteration.txt")

  create_file_if(c(node_it, highest_it_file), keep = TRUE)
  create_file_if(c(ram_usage, err_file, warn_file))
  
  if (!is.null(iterations)) {
    
    batch_iterations <- iterations
    
    catn("Initiating specific iteration(s)", highcat(batch_iterations))
    
  } else {
    node_its <- readLines(node_it_file)
    vebcat("Node iterations:", node_its, veb = verbose)
    
    if (!is.na(node_its[1])) {
      node_it <- gsub("node", "", node_its)
      
      vebcat("Node iterations from previous session:", node_it, veb = verbose)
      
      start_iteration <- as.integer(min(node_it)) # This one will not have run
      
      catn("Start iteration:", highcat(start_iteration))
    } else {
      catn("Node iterations file is null.")
      
      highest_it <- as.integer(readLines(highest_it_file))
      if (is.na(highest_it[1])) {
        catn("Highest iteration is null. Assuming sequence has never been run before.")
        start_iteration <- 0
      } else {
        catn("Previous session completed successfully on iteration", highcat(highest_it))
        catn("Input list is expected to take", highcat(length(spec.list)), "iterations.")
        start_iteration <- highest_it + 1
      }
    }
    
    if (start_iteration >= length(spec.list)) {
      catn("Start iteration:", highcat(start_iteration), "number of species:", highcat(length(spec.list)))
      stop("STOP: Previous iteration is higher or the same as the number of species.")
    }
    
    # If iterations is not provided, start from the highest saved iteration
    i <- start_iteration
    end <- length(spec.list)
    batch_iterations <- i:end
    
    catn("Initiating from iteration:", highcat(i), "to", highcat(end))
  }
  
  cores_max <- min(length(batch_iterations), cores.max)
  
  catn("Creating cluster of", cores_max, "cores.", veb = verbose)
  
  cl <- makeCluster(cores.max)
  
  vebcat("Including the necessary components in each core.", veb = verbose)
  
  cluster_params <- c(
    "par.dir", 
    "spec.list", 
    "iterations", 
    "cores.max", 
    "cores.max.high", 
    "min.disk.space", 
    "verbose",
    names(custom.exports)
  )
  
  clusterExport(cl, cluster_params, envir = environment())
  
  clusterEvalQ(cl, {
    for (file in cluster.evals) {
      source(file)
    }
  })
  
  vebcat("Creating a vector for the results.", veb = verbose)
  
  current_disk_space <- get_disk_space("/export", units = "GB")
  
  catn("\nRemaining disk space (GB)")
  cat(sprintf("%8s | %8s | %8s \n", "Minimum", "Current", "Remaining"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", min.disk.space, current_disk_space, current_disk_space - min.disk.space))
  
  catn("\nMemory allocation (GB)")
  cat(sprintf("%8s | %8s | %8s \n", "Maximum", "Limit", "Current"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", mem_total / 1024^3, mem_limit / 1024^3, get_mem_usage("used", format = "gb")))
  
  catn("Hypervolume sequence has started, progress is being logged to:", colcat(logs_dir, color = "output"))

  return(list(
    ram.use = ram_usage,
    batch = batch_iterations,
    highest.iteration = highest_it_file
  ))
}
