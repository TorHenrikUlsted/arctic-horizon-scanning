parallell_processing <- function(sp_list, method, accuracy, project, proj.incl.t, iterations = NULL, max_cores, min.disk.space, show.plot, verbose) {
  on.exit(closeAllConnections())
  
  parallell_timer <- start_timer("parallell_timer")

  cat(blue("Initiating hypervolume sequence \n"))
  cat(cc$aquamarine("Outputs are being logged at ouputs/hypervolume/logs. \n"))

  if (!file.exists(paste0("./outputs/hypervolume/logs/", method, "-error.txt"))) file.create(paste0("./outputs/hypervolume/logs/", method, "-error.txt"))
  if (!file.exists(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"))) file.create(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"))

  stats_dir <- paste0("./outputs/hypervolume/stats")
  create_dir_if(stats_dir)

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

  cl <- makeCluster(max_cores)

  highest_it_file <- "./outputs/hypervolume/logs/highest_iteration.txt"

  if (is.null(iterations)) {
    if (file.exists(highest_it_file)) {
      highest_iteration <- as.integer(readLines(highest_it_file))

      cat("Previous session found, continuing from iteration:", cc$lightSteelBlue(highest_iteration + 1), "with", cc$lightSteelBlue(max_cores), "cores. \n")

      if (highest_iteration >= length(sp_list)) {
        stop(cc$lightCoral("STOP: Previous iteration is higher or the same as the number of species. \n"))
      }
    } else {
      highest_iteration <- 0
    }
    
    # If iterations is not provided, start from the highest saved iteration
    i <- highest_iteration + 1
    end <- length(sp_list)
    batch_iterations <- i:end
  } else {
    batch_iterations <- iterations

    cat("Initiating", cc$lightSteelBlue(length(batch_iterations)), "specific iteration(s) with", cc$lightSteelBlue(max_cores), "cores. \n")
  }

  clusterEvalQ(cl, {
    source("./src/utils/utils.R")
    source("./src/setup/setup.R")
    source("./src/hypervolume/hypervolume.R")
  })

  results <- vector("list", length(sp_list))

  current_disk_space <- get_disk_space("/home", units = "GB")

  cat("Remaining disk space (GB) \n")
  cat(sprintf("%8s | %8s | %8s \n", "Minimum", "Current", "Remaining"))
  cat(sprintf("%8.2f | %8.2f | %8.0f \n", min.disk.space, current_disk_space, current_disk_space - min.disk.space))

  res <- clusterApplyLB(cl, batch_iterations, function(j) {
    node_processing(j, sp_list, proj.incl.t, method, accuracy, project, show.plot, verbose, min.disk.space)
  })

  results[batch_iterations] <- res


  cat("Finishing up \n")

  stopCluster(cl)

  highest_iteration <- max(unlist(lapply(results, function(res) res$iteration)))

  writeLines(as.character(highest_iteration), highest_it_file)
  
  end_timer(parallell_timer)

  cat(cc$lightGreen("Hypervolume sequence completed succesfully \n"))
}
