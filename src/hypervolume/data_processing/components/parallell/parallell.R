parallell_processing <- function(spec.list, method, accuracy, hv.projection, proj.incl.t, iterations = NULL, cores.max.high = 1, cores.max = 1, min.disk.space, hv.dir, show.plot = F, verbose = T) {
  on.exit(closeAllConnections())
  
  cat(blue("Initiating hypervolume sequence \n"))
  
  parallell_timer <- start_timer("parallell_timer")

  if (verbose) cat("Setting up hypervolume files and folders. \n")
  
  logs_dir <-   paste0(hv.dir, "/logs")
  proj_dir <-   paste0(hv.dir, "/projections")
  stats_dir <- paste0(hv.dir, "/stats")
  locks_dir <- paste0(hv.dir, "/locks")
  
  create_dir_if(logs_dir)
  create_dir_if(proj_dir)
  create_dir_if(stats_dir)
  create_dir_if(locks_dir)
  
  highest_it_file <- paste0(hv.dir, "/logs/highest_iteration.txt")
  ram_usage <- paste0(logs_dir, "/", method, "-ram-usage.txt")
  err_file <- paste0(logs_dir, "/", method, "-error.txt")
  warn_file <- paste0(logs_dir, "/", method, "-warning.txt")
  node_it <- paste0(logs_dir, "/node-iterations.txt")
  
  create_file_if(node_it, keep = TRUE)
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
    if (file.exists(highest_it_file)) {
      highest_iteration <- as.integer(readLines(highest_it_file))

      cat("Previous session found, continuing from iteration:", cc$lightSteelBlue(highest_iteration + 1), "with", cc$lightSteelBlue(cores.max), "core(s). \n")

      if (highest_iteration >= length(spec.list)) {
        stop(cc$lightCoral("STOP: Previous iteration is higher or the same as the number of species. \n"))
      }
    } else {
      create_file_if(highest_it_file, keep = TRUE)
      highest_iteration <- 0
    }
    
    # If iterations is not provided, start from the highest saved iteration
    i <- highest_iteration + 1
    end <- length(spec.list)
    batch_iterations <- i:end
  } else {
    batch_iterations <- iterations

    cat("Initiating", cc$lightSteelBlue(length(batch_iterations)), "specific iteration(s) with", cc$lightSteelBlue(cores.max), "core(s),", "and a high max of", cc$lightSteelBlue(cores.max.high), "\n")
  }
  
  if (verbose) cat("Creating cluster of", cc$lightSteelBlue(cores.max), "core(s), with a high max of", cc$lightSteelBlue(cores.max.high), "\n")
  
  cl <- makeCluster(cores.max)
  
  clusterExport(cl, c("spec.list", "proj.incl.t", "method", "accuracy", "hv.projection", "cores.max.high", "min.disk.space", "hv.dir", "show.plot", "verbose"), envir = environment())
  
  clusterEvalQ(cl, {
    source("./src/utils/utils.R")
    source("./src/setup/setup.R")
    source("./src/hypervolume/hypervolume.R")
  })
  
  if (verbose) cat("Creating a vector for the results. \n")

  results <- vector("list", length(spec.list))

  current_disk_space <- get_disk_space("/export", units = "GB")
  
  cat("\nRemaining disk space (GB) \n")
  cat(sprintf("%8s | %8s | %8s \n", "Minimum", "Current", "Remaining"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", min.disk.space, current_disk_space, current_disk_space - min.disk.space))
  
  cat("\nMemory allocation (GB) \n")
  cat(sprintf("%8s | %8s | %8s \n", "Maximum", "Current", "Limit"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", mem_total / 1024^3, get_mem_usage("free", format = "gb"), mem_limit / 1024^3))
  
  cat("Hypervolume sequence has started, progress is being logged to:", yellow(logs_dir), "\n")
  
  res <- clusterApplyLB(cl, batch_iterations, function(j) {
    ram_msg <- FALSE
    # RAM check
    mem_used_gb <- get_mem_usage(type = "used", format = "gb")
    mem_limit_gb <- mem_limit / 1024^3
    
    while (mem_used_gb >= mem_limit_gb) {
      if (!ram_msg) {
        ram_con <- file(ram_usage, open = "a")
        writeLines(paste0("RAM usage",mem_used_gb,"is above the maximum",mem_limit_gb,"Waiting with node", j), ram_con)
        close(ram_con)
      }
      Sys.sleep(5)  # Wait for 5 seconds before checking again
      Sys.sleep(runif(1, 0, 1)) # Add random seconds between 0 and 1 to apply difference if multiple nodes are in queue
      mem_used_gb <- get_mem_usage(type = "used", format = "gb")
    }

    
    node_processing(j, spec.list, proj.incl.t, method, accuracy, hv.projection, cores.max.high, min.disk.space, hv.dir, show.plot, verbose)
    
  })

  results[batch_iterations] <- res

  cat("Finishing up \n")

  stopCluster(cl)

  highest_iteration <- max(unlist(lapply(results, function(res) res$iteration)))

  writeLines(as.character(highest_iteration), highest_it_file)
  
  end_timer(parallell_timer)

  cat(cc$lightGreen("Hypervolume sequence completed succesfully \n"))
}
