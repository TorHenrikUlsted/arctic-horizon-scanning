hypervolume_sequence <- function(spec.list, method, accuracy, hv.dims, hv.projection, proj.incl.t, iterations = NULL, cores.max.high = 1, cores.max = 1, min.disk.space, hv.dir, show.plot = F, verbose = T) {
  on.exit(closeAllConnections())
  
  vebcat("Initiating hypervolume sequence", color = "seqInit")
  
  parallell_timer <- start_timer("parallell_timer")
  
  vebcat("Setting up hypervolume files and folders.", veb = verbose)
  
  logs_dir <-   paste0(hv.dir, "/logs")
  proj_dir <-   paste0(hv.dir, "/projections")
  stats_dir <- paste0(hv.dir, "/stats")
  locks_dir <- paste0(hv.dir, "/locks")
  
  create_dir_if(logs_dir)
  create_dir_if(proj_dir)
  create_dir_if(stats_dir)
  create_dir_if(locks_dir)
  
  ram_usage <- paste0(logs_dir, "/", method, "-ram-usage.txt")
  err_file <- paste0(logs_dir, "/", method, "-error.txt")
  warn_file <- paste0(logs_dir, "/", method, "-warning.txt")
  node_it <- paste0(logs_dir, "/node-iterations.txt")
  highest_it_file <- paste0(logs_dir, "/highest_iteration.txt")
  
  create_file_if(node_it, keep = TRUE)
  create_file_if(highest_it_file, keep = TRUE)
  create_file_if(ram_usage)
  create_file_if(err_file)
  create_file_if(warn_file)

  if (!file.exists(paste0(stats_dir, "/", method, "-stats.csv"))) {
    df <- data.frame(
      species = character(0),
      observations = integer(0),
      dimensions = integer(0),
      samplesPerPoint = integer(0),
      randomPoints = integer(0),
      excluded = logical(0),
      jaccard = numeric(0),
      sorensen = numeric(0),
      fracVolumeSpecies = numeric(0),
      fracVolumeRegion = numeric(0),
      overlapRegion = numeric(0),
      includedOverlap = numeric(0)
    )

    fwrite(df, paste0(stats_dir, "/", method, "-stats.csv"), row.names = F, bom = T)
  }

  if (is.null(iterations)) {
      node_its <- readLines(node_it)
      catn("Node iterations:", node_its)
      
      if (is.na(node_its[1])) {
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
        
      } else {
        node_int <- gsub("node", "", node_its)
        
        catn("Node iterations from previous session:", node_int)
        
        start_iteration <- as.integer(min(node_int))
        
        catn("Start iteration:", highcat(start_iteration))
      }
      
      if (start_iteration >= length(spec.list)) {
        catn("Start iteration:", highcat(start_iteration), "number of species:", highcat(length(spec.list)))
        stop("STOP: Previous iteration is higher or the same as the number of species.")
      }
      
      # If iterations is not provided, start from the highest saved iteration
      i <- start_iteration
      end <- length(spec.list)
      batch_iterations <- i:end
      
      cores.max <- min(length(batch_iterations), cores.max)
      cores.max.high <- min(cores.max, cores.max.high)
      
      catn("Initiate from iteration:", highcat(i), "with", highcat(cores.max), "core(s), and a high max of", highcat(cores.max.high), "\n")
    
  } else {
    batch_iterations <- iterations
    
    cores.max <- min(length(batch_iterations), cores.max)
    cores.max.high <- min(cores.max, cores.max.high)

    catn(
      "Initiating specific iteration(s)", highcat(batch_iterations), 
      "with", highcat(cores.max), "core(s),", 
      "and a high max of", highcat(cores.max.high)
    )
    
  }
  
  vebcat("Creating cluster of cores.", veb = verbose)
  
  cl <- makeCluster(cores.max)
  
  vebcat("Including the necessary components in each core.", veb = verbose)
  
  clusterExport(cl, c("spec.list", "proj.incl.t", "method", "accuracy", "hv.dims", "hv.projection", "cores.max.high", "min.disk.space", "hv.dir", "show.plot", "verbose"), envir = environment())
  
  clusterEvalQ(cl, {
    source("./src/utils/utils.R")
    source("./src/setup/setup.R")
    source("./src/hypervolume/data_acquisition/data_acquisition.R")
    source("./src/hypervolume/data_analysis/data_analysis.R")
    source("./src/hypervolume/data_processing/data_processing.R")
    source("./src/hypervolume/hv_analysis/hv_analysis.R")
    source("./src/hypervolume/node.R")
  })
  
  vebcat("Creating a vector for the results.", veb = verbose)

  results <- vector("list", length(spec.list))

  current_disk_space <- get_disk_space("/export", units = "GB")
  
  catn("\nRemaining disk space (GB)")
  cat(sprintf("%8s | %8s | %8s \n", "Minimum", "Current", "Remaining"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", min.disk.space, current_disk_space, current_disk_space - min.disk.space))
  
  catn("\nMemory allocation (GB)")
  cat(sprintf("%8s | %8s | %8s \n", "Maximum", "Limit", "Current"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", mem_total / 1024^3, mem_limit / 1024^3, get_mem_usage("used", format = "gb")))
  
  catn("Hypervolume sequence has started, progress is being logged to:", colcat(logs_dir, color = "output"))
  
  res <- clusterApplyLB(cl, batch_iterations, function(j) {
    ram_msg <- FALSE
    # RAM check
    mem_used_gb <- get_mem_usage(type = "used", format = "gb")
    mem_limit_gb <- mem_limit / 1024^3
    
    while (mem_used_gb >= mem_limit_gb) {
      if (!ram_msg) {
        ram_con <- file(ram_usage, open = "a")
        writeLines(paste0("RAM usage ", mem_used_gb, " is above the maximum ", mem_limit_gb, " Waiting with node", j), ram_con)
        close(ram_con)
        ram_msg = TRUE
      }
      Sys.sleep(60)  # Wait for 5 seconds before checking again
      Sys.sleep(runif(1, 0, 1)) # Add random seconds between 0 and 1 to apply difference if multiple nodes are in queue
      mem_used_gb <- get_mem_usage(type = "used", format = "gb")
    }
    
    node_processing(j, spec.list, proj.incl.t, method, accuracy, hv.dims, hv.projection, cores.max.high, min.disk.space, hv.dir, show.plot, verbose)
    
  })

  results[batch_iterations] <- res

  catn("Finishing up.")

  stopCluster(cl)
  
  prev_highest_it <- as.integer(readLines(highest_it_file))

  highest_iteration <- max(unlist(lapply(results, function(res) res$iteration)))
  
  if (highest_iteration > prev_highest_it) writeLines(as.character(highest_iteration), highest_it_file)
  
  end_timer(parallell_timer)

  veb("Hypervolume sequence completed succesfully", color = "seqSuccess")
}
