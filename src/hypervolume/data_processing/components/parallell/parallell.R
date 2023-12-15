parallell_processing <- function(sp_list, method, max_cores, projection, show_plot, verbose, iterations = NULL) {
  on.exit(closeAllConnections())
  
  cat(blue("Initiating hypervolume sequence \n"))
  cat(cc$aquamarine("Outputs are being logged at ouputs/hypervolume/logs. \n"))

  highest_it_file <- "./outputs/hypervolume/logs/highest_iteration.txt"
  
  if (!file.exists(paste0("./outputs/hypervolume/logs/", method, "-error.txt"))) file.create(paste0("./outputs/hypervolume/logs/", method, "-error.txt"))
  if (!file.exists(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"))) file.create(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"))
  stats_dir <- paste0("./outputs/hypervolume/stats")
  create_dir_if(stats_dir)
  if (!file.exists(paste0(stats_dir, "/", method, "-stats.csv"))) {
    df <- data.frame(
      Species = character(0),
      jaccard = numeric(0),
      sorensen = numeric(0),
      fracVolumeSpecies = numeric(0),
      fracVolumeRegion = numeric(0),
      overlapRegion = numeric(0),
      randomPoints = integer(0),
      samplesPerPoint = integer(0),
      observations = integer(0),
      dimensions = integer(0),
      maxInclusionPoints = integer(0),
      includedPoints = integer(0),
      excludedPoints = integer(0),
      includedOverlap = numeric(0)
    )
  fwrite(df, paste0(stats_dir, "/", method, "-stats.csv"), row.names = F, bom = T)
  }
  
  min_disk_space <- get_disk_space("/home", units = "GB") * 0.2

  cl <- makeCluster(max_cores)
  
  if (file.exists(highest_it_file)) {
    highest_iteration <- as.integer(readLines(highest_it_file))
    cat("Previous session found, continuing from iteration:", cc$lightSteelBlue(highest_iteration), "with", cc$lightSteelBlue(max_cores), "cores. \n")
    
    if (highest_iteration >= length(sp_list)) {
      stop(cc$lightCoral("STOP: Previous iteration is higher or the same as the length of species. \n"))
    }
  } else {
    highest_iteration <- 0
  }
  
  if (is.null(iterations)) {
    # If iterations is not provided, start from the highest saved iteration
    i <- highest_iteration + 1
    end <- length(sp_list)
    batch_iterations <- i:end
  } else {
    # If iterations is provided, use it
    i <- min(iterations)
    end <- max(iterations)
    batch_iterations <- i:end
  }

  clusterEvalQ(cl, {
    source("./src/utils/utils.R")
    source("./src/setup/setup.R")
    source("./src/hypervolume/hypervolume.R")
  })

  results <- vector("list", length(sp_list))
  
  while (i <= end) {
    current_disk_space <- get_disk_space("/home", units = "GB")
    
    cat("Remaining disk space \n")
    cat(sprintf("%8s | %8s | %8s \n", "Minimum", "Current", "Percentage"))
    cat(sprintf("%8.4f | %8.4f | %8.0f %% \n", min_disk_space, current_disk_space, (current_disk_space - min_disk_space) / current_disk_space * 100 ))

    if (current_disk_space <= min_disk_space) {
      cat("Disk space has reached 80% of the total. No new iterations will be started.\n")
      break
    }
      
      res <- clusterApplyLB(cl, batch_iterations, function(j) {
        node_processing(j, sp_list, method, show_plot, verbose)
      })

      results[batch_iterations] <- res


   i <- i + length(batch_iterations)
    
  }
  
  cat("Finishing up \n")
  
  stopCluster(cl)

  highest_iteration <- max(unlist(lapply(results, function(res) res$iteration)))

  writeLines(as.character(highest_iteration), highest_it_file)

  cat(cc$lightGreen("Hypervolume sequence completed succesfully \n"))
}
