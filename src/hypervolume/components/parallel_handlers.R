setup_parallel <- function(par.dir, spec.list, iterations, cores.max, cores.max.high, min.disk.space, verbose = FALSE, custom.exports, custom.evals) {
  vebcat("Setting up hypervolume files and folders.", veb = verbose)

  logs_dir <- paste0(par.dir, "/logs")

  create_dir_if(logs_dir)

  ram_usage <- paste0(logs_dir, "/ram-usage.txt")
  node_it_file <- paste0(logs_dir, "/node-iterations.txt")
  highest_it_file <- paste0(logs_dir, "/highest-iteration.txt")
  warn_file <- paste0(logs_dir, "/warning.txt")
  err_file <- paste0(logs_dir, "/error.txt")

  create_file_if(node_it_file, highest_it_file, keep = TRUE)
  create_file_if(ram_usage, warn_file, err_file)
  finished <- FALSE

  if (!is.null(iterations)) {
    batch_iterations <- iterations

    catn("Initiating specific iteration(s)", highcat(batch_iterations))
  } else {
    node_its <- readLines(node_it_file)
    vebcat("Node iterations:", node_its, veb = verbose)

    if (is.na(node_its[1]) || node_its[1] == "") {
      catn("Node iterations file is empty.")

      highest_it <- as.integer(readLines(highest_it_file))
      if (is.na(highest_it[1])) {
        catn("Highest iteration is null. Assuming sequence has never been run before.")
        start_iteration <- 1
      } else {
        catn("Previous session completed successfully on iteration", highcat(highest_it))
        catn("Input list is expected to take", highcat(length(spec.list)), "iterations.")
        start_iteration <- highest_it + 1
      }
    } else {
      node_it <- gsub("node", "", node_its)

      vebcat("Node iterations from previous session:", node_it, veb = verbose)

      start_iteration <- as.integer(min(node_it))

      catn("Previous session found, Continuing from iteration", highcat(start_iteration))
    }

    if (start_iteration >= length(spec.list)) {
      catn("Start iteration:", highcat(start_iteration), "number of species:", highcat(length(spec.list)))
      finished <- TRUE
      return(list(
        finished = finished
      ))
    }

    # If iterations is not provided, start from the highest saved iteration
    i <- start_iteration
    end <- length(spec.list)
    batch_iterations <- i:end

    catn("Initiating from iteration:", highcat(i), "to", highcat(end))
  }

  cores_max <- min(length(batch_iterations), cores.max)

  catn("Creating cluster of", highcat(cores_max), "core(s).")

  cl <- makeCluster(cores_max)

  vebcat("Including the necessary components in each core.", veb = verbose)

  cluster_params <- c(
    "par.dir",
    "spec.list",
    "iterations",
    "cores.max",
    "cores.max.high",
    "min.disk.space",
    "verbose",
    names(custom.exports),
    "custom.evals"
  )

  clusterExport(cl, cluster_params, envir = environment())

  clusterEvalQ(cl, {
    for (file in custom.evals) {
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
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", config$memory$mem_total / 1024^3, config$memory$config$memory$mem_limit / 1024^3, get_mem_usage("used", format = "gb")))

  catn("Hypervolume sequence has started, progress is being logged to:", colcat(logs_dir, color = "output"))

  return(list(
    cl = cl,
    cores = cores_max,
    ram.use = ram_usage,
    batch = batch_iterations,
    highest.iteration = highest_it_file,
    finished = finished
  ))
}

optimize_queue <- function(dt, cores.max, high_ram_threshold = 0.2, verbose = FALSE) {
  setorder(dt, -medianLat)
  n_high_ram <- ceiling(nrow(dt) * high_ram_threshold)
  high_ram_species <- dt[1:n_high_ram, filename]
  low_ram_species <- dt[(n_high_ram + 1):nrow(dt), filename]
  
  chunks <- vector("list", cores.max)
  
  for (i in seq_along(high_ram_species)) {
    core_index <- (i - 1) %% total_cores + 1
    chunks[[core_index]] <- c(chunks[[core_index]], high_ram_species[i])
  }
  
  current_core <- 1
  for (species in low_ram_species) {
    while (length(chunks[[current_core]]) >= ceiling(length(high_ram_species) / total_cores)) {
      current_core <- (current_core %% total_cores) + 1
    }
    chunks[[current_core]] <- c(chunks[[current_core]], species)
    current_core <- (current_core %% total_cores) + 1
  }
  
  # Print statistics
  vebcat("Created", highcat(length(chunks)), "chunks")
  vebcat("High-RAM species:", highcat(length(high_ram_species)))
  vebcat("Low-RAM species:", highcat(length(low_ram_species)))
  
  # Print distribution of high-RAM species across chunks
  high_ram_dist <- sapply(chunks, function(chunk) sum(chunk %in% high_ram_species))
  vebprint(high_ram_dist, text = "Distribution of high-RAM species across chunks:")
  
  # Print average median latitude for each chunk
  avg_lat <- sapply(chunks, function(chunk) {
    mean(spec_count_dt[filename %in% chunk, medianLat])
  })
  
  vebcat("Average median latitude for each chunk:", highcat(round(avg_lat, 2)))
  return(chunks) 
}
