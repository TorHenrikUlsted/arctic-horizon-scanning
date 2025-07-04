setup_parallel <- function(par.dir, spec.list, iterations, cores.max, cores.max.high, min.disk.space, verbose = FALSE, custom.exports, custom.evals) {
  vebcat("Setting up hypervolume files and folders.", veb = verbose)

  logs_dir <- paste0(par.dir, "/logs")

  create_dir_if(logs_dir)

  ram_usage <- paste0(logs_dir, "/ram-usage.txt")
  node_it_file <- paste0(logs_dir, "/node-iterations.txt")
  highest_it_file <- paste0(logs_dir, "/highest-iteration.txt")
  warn_file <- paste0(logs_dir, "/warning.txt")
  err_file <- paste0(logs_dir, "/error.txt")
  time_est <- paste0()

  create_file_if(node_it_file, highest_it_file, warn_file, err_file, keep = TRUE)
  create_file_if(ram_usage, keep = FALSE)
  finished <- FALSE

  if (!is.null(iterations)) {
    batch_iterations <- iterations
    catn("Initiating specific iteration(s)", highcat(batch_iterations))
  } else {
    node_its <- as.integer(gsub("node", "", readLines(node_it_file)))
    highest_it <- as.integer(readLines(highest_it_file))

    vebcat("Node iterations:", node_its, veb = verbose)
    vebcat("Highest completed iteration:", highest_it, veb = verbose)

    if (length(node_its) == 0 || is.na(node_its[1])) {
      start_iteration <- if (length(highest_it) == 0) 1 else highest_it + 1
    } else {
      start_iteration <- max(highest_it, max(node_its)) + 1
    }

    if (start_iteration > length(spec.list) & length(node_its) == 0) {
      catn("All iterations completed. No more processing needed.")
      return(list(finished = TRUE))
    }

    end <- length(spec.list)

    batch_iterations <- c(
      if (length(node_its) > 0) unique(node_its),
      if (start_iteration <= end) start_iteration:end
    )

    if (start_iteration > length(spec.list)) {
      catn("Iterations finished")
    } else {
      catn("Processing from iteration:", highcat(start_iteration), "/", highcat(end))
    }

    if (length(node_its) != 0) catn("Including node iterations:", highcat(paste(unique(node_its), collapse = ", ")))
  }

  cores_max <- min(length(batch_iterations), cores.max)

  catn("Creating cluster of", highcat(cores_max), "core(s).")

  cl <- makeCluster(cores_max)

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

  vebcat("Exporting cluster parameters", veb = verbose)

  clusterExport(cl, cluster_params, envir = environment())

  vebcat("Including the necessary components in each core.", veb = verbose)

  tryCatch(
    {
      clusterEvalQ(cl, {
        source("./R/utils/utils.R")
        load_utils(parallel = TRUE)
        for (file in custom.evals) {
          tryCatch(
            {
              source(file, local = TRUE) # Use local=TRUE to avoid polluting global env
            },
            error = function(e) {
              stop(paste("Error in file", file, ":", e$message))
            }
          )
        }
        gc(full = TRUE)
      })
    },
    error = function(e) {
      vebcat("Error when sourcing files for each core.", color = "fatalError")
      vebprint(custom.evals, text = "Files attempted to source:")
      stopCluster(cl)
      closeAllConnections()
      stop(e)
    }
  )

  vebcat("Creating a vector for the results.", veb = verbose)

  current_disk_space <- get_disk_space("/export", units = "GB")

  catn("\nRemaining disk space (GB)")
  cat(sprintf("%8s | %8s | %8s \n", "Minimum", "Current", "Remaining"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", min.disk.space, current_disk_space, current_disk_space - min.disk.space))

  catn("\nMemory allocation (GB)")
  cat(sprintf("%8s | %8s | %8s \n", "Maximum", "Limit", "Current"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", config$memory$mem_total / 1024^3, config$memory$mem_limit / 1024^3, get_mem_usage("used", format = "gb")))

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
  catn("Optimizing queue")

  vebcat("Input data:", veb = verbose)
  vebprint(dt, veb = verbose)

  if (nrow(dt) <= cores.max) {
    vebcat("Number of species is less than or equal to number of cores. Returning one species per chunk.", veb = verbose)
    chunks <- lapply(1:nrow(dt), function(i) dt$filename[i])
    return(chunks)
  }

  setorder(dt, -medianLat)
  n_high_ram <- max(1, min(ceiling(nrow(dt) * high_ram_threshold), nrow(dt) - 1))
  high_ram_species <- dt[1:n_high_ram, filename]
  low_ram_species <- dt[(n_high_ram + 1):nrow(dt), filename]

  vebcat("High RAM species:", veb = verbose)
  vebprint(high_ram_species, veb = verbose)
  vebcat("Low RAM species:", veb = verbose)
  vebprint(low_ram_species, veb = verbose)

  chunks <- vector("list", cores.max)

  for (i in seq_along(high_ram_species)) {
    core_index <- (i - 1) %% cores.max + 1
    chunks[[core_index]] <- c(chunks[[core_index]], high_ram_species[i])
  }

  vebcat("Chunks after distributing high RAM species:", veb = verbose)
  vebprint(chunks, veb = verbose)

  if (length(low_ram_species) > 0) {
    current_core <- 1
    for (species in low_ram_species) {
      chunks[[current_core]] <- c(chunks[[current_core]], species)
      current_core <- (current_core %% cores.max) + 1
    }
  }

  # Remove empty chunks
  chunks <- chunks[sapply(chunks, length) > 0]

  # Print statistics
  vebcat("Created", highcat(length(chunks)), "chunks", veb = verbose)
  vebcat("High-RAM species:", highcat(length(high_ram_species)), veb = verbose)
  vebcat("Low-RAM species:", highcat(length(low_ram_species)), veb = verbose)

  # Print distribution of high-RAM species across chunks
  high_ram_dist <- sapply(chunks, function(chunk) sum(chunk %in% high_ram_species))
  vebprint(high_ram_dist, text = "Distribution of high-RAM species across chunks:")

  # Print average median latitude for each chunk
  avg_lat <- sapply(chunks, function(chunk) {
    mean(dt[filename %in% chunk, medianLat])
  })

  vebprint(round(avg_lat, 2), text = "Average median latitude for each chunk:")
  return(chunks)
}

#' Clean up files in a directory
#'
#' @param dir_path Character string specifying the directory path
#' @param recursive Logical, whether to clean files in subdirectories (default: FALSE)
#' @param remove_dirs Logical, whether to remove directories too (default: FALSE)
#' @param pattern Character string containing a regular expression to match files (default: NULL)
#' @param verbose Logical, whether to print progress messages (default: FALSE)
#' @return Named logical vector indicating which files were successfully removed
#' @examples
#' \dontrun{
#' # Clean only files in current directory
#' cleanup_files("./temp")
#'
#' # Clean files including subdirectories
#' cleanup_files("./temp", recursive = TRUE)
#'
#' # Remove everything including directories
#' cleanup_files("./temp", recursive = TRUE, remove_dirs = TRUE)
#'
#' # Clean only specific files
#' cleanup_files("./temp", pattern = "\\.txt$")
#' }
cleanup_files <- function(dir_path,
                          recursive = FALSE,
                          remove_dirs = FALSE,
                          pattern = NULL,
                          verbose = FALSE) {
  # Check if directory exists
  if (!dir.exists(dir_path)) {
    warning(sprintf("Directory does not exist: %s", dir_path))
    return(invisible(logical(0)))
  }

  # List all files
  files <- list.files(dir_path,
    recursive = recursive,
    full.names = TRUE,
    all.files = TRUE, # Include hidden files
    include.dirs = remove_dirs,
    pattern = pattern
  )

  if (length(files) == 0) {
    if (verbose) catn("No files found to remove.")
    return(invisible(logical(0)))
  }

  # If not removing directories, filter them out
  if (!remove_dirs) {
    files <- files[!dir.exists(files)]
  } else {
    # If removing directories, sort them to remove deepest first
    files <- sort(files, decreasing = TRUE)
  }

  if (verbose) {
    catn(sprintf("Found %d items to remove in %s", length(files), dir_path))
    if (length(files) > 1) {
      catn("Files to remove:")
      print(files)
    }
  }

  # Try to remove each file
  removed <- logical(length(files))
  names(removed) <- files

  for (i in seq_along(files)) {
    if (remove_dirs && dir.exists(files[i])) {
      # Use unlink for directories
      removed[i] <- !as.logical(unlink(files[i], recursive = TRUE))
    } else {
      # Use file.remove for files
      removed[i] <- file.remove(files[i])
    }

    if (verbose && !removed[i]) {
      warning(sprintf("Failed to remove: %s", files[i]))
    }
  }

  if (verbose) {
    catn(sprintf(
      "Successfully removed %d/%d items",
      sum(removed),
      length(removed)
    ))
  }

  invisible(removed)
}
