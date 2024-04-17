hypervolume_sequence <- function(
    spec.list,
    iterations = NULL, 
    cores.max.high = 1,
    cores.max = 1,
    min.disk.space = 1,
    verbose = FALSE,
    hv.dir,
    hv.method, 
    hv.accuracy, 
    hv.dims, 
    hv.incl.threshold
    ) {
  on.exit(closeAllConnections())
  
  vebcat("Initiating hypervolume sequence", color = "seqInit")
  
  parallell_timer <- start_timer("parallell_timer")
  
  hv_dir <- paste0(hv.dir, "/", hv.method, "-sequence")
  
  create_dir_if(hv_dir)
  
  stats_file <- paste0(hv_dir, "/stats.csv")
  removed_file <- paste0(hv_dir, "/removed-species.csv")
  
  init_dt <- data.table(
    cleanName = character(0),
    iteration = integer(0),
    observations = integer(0),
    dimensions = integer(0),
    samplesPerPoint = integer(0),
    randomPoints = integer(0),
    excluded = logical(0),
    jaccard = numeric(0),
    sorensen = numeric(0),
    fracVolumeSpecies = numeric(0),
    fracVolumeRegion = numeric(0),
    realizedNiche = numeric(0),
    overlapRegion = numeric(0),
    includedOverlap = numeric(0),
    kingdom = character(0), 
    phylum = character(0), 
    class = character(0),
    order = character(0),
    family = character(0), 
    genus = character(0), 
    species = character(0),
    infraspecificEpithet = character(0),
    taxonRank = character(0), 
    scientificName = character(0),
    countryCode = character(0),
    country = character(0),
    meanLong = numeric(0),
    meanLat = numeric(0)
  )
  
  if (!file.exists(stats_file)) {
    fwrite(init_dt, stats_file, row.names = F, bom = T)
  }
  
  custom_exports = c(
    hv.dir,
    hv.method, 
    hv.accuracy, 
    hv.dims, 
    hv.incl.threshold
  )
  
  custom_evals = c(
    "./src/utils/utils.R",
    "./src/setup/setup_sequence.R",
    "./src/hypervolume/hypervolume.R",
    "./src/hypervolume/node_hypervolume.R",
    "./src/hypervolume/parallel_hypervolume.R"
  )
  
  spec_count_dt <- count_observations(
    spec.list = spec.list,
    dimensions = hv.dims,
    verbose = verbose
  )
  
  spec_removed <- spec_count_dt[removed == TRUE, ]
  spec_removed[, filename := NULL]

  fwrite(spec_removed, removed_file, bom = TRUE)

  spec_count_dt <- spec_count_dt[removed == FALSE, ]

  setorder(spec_count_dt, meanLat)
  
  print(spec_count_dt)

  spec_list <- spec_count_dt$filename

  vebprint(head(spec_list, 10), verbose, "species list sample:")
  
  stop()
  
  parallel <- setup_parallel(
    par.dir = hv_dir,
    spec.list = spec.list,
    iterations = iterations,
    cores.max = cores.max,
    cores.max.high = cores.max.high,
    min.disk.space = min.disk.space,
    verbose = verbose,
    custom.exports = custom_exports,
    custom.evals = custom_evals
  )
  
  results <- vector("list", length(spec.list))
  
  vebprint(clusterEvalQ(parallel$cl, ls()), veb = verbose, text = "All cluster variables:")

  res <- clusterApplyLB(parallel$cl, parallel$batch, function(j) {
    ram_msg <- FALSE
    # RAM check
    mem_used_gb <- get_mem_usage(type = "used", format = "gb")
    mem_limit_gb <- mem_limit / 1024^3
    
    while (mem_used_gb >= mem_limit_gb) {
      if (!ram_msg) {
        ram_con <- file(parallel$ram.use, open = "a")
        writeLines(paste0("RAM usage ", mem_used_gb, " is above the maximum ", mem_limit_gb, " Waiting with node", j), ram_con)
        close(ram_con)
        ram_msg = TRUE
      }
      Sys.sleep(60)  # Wait for 5 seconds before checking again
      Sys.sleep(runif(1, 0, 1)) # Add random seconds between 0 and 1 to apply difference if multiple nodes are in queue
      mem_used_gb <- get_mem_usage(type = "used", format = "gb")
    }
    
    node_hypervolume(
      process.dir = par.dir,
      iteration = j,
      spec.list = spec.list,
      columns.to.read = c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "infraspecificEpithet",
        "taxonRank",
        "scientificName",
        "cleanName",
        "decimalLongitude", 
        "decimalLatitude", 
        "coordinateUncertaintyInMeters", 
        "countryCode",
        "occurrenceStatus",
        "stateProvince",
        "year"
      ),
      min.disk.space = min.disk.space,
      cores.max.high = cores.max.high,
      init.dt = init_dt,
      verbose = verbose,
      hv.incl.threshold = hv.incl.threshold, 
      hv.method = hv.method, 
      hv.accuracy = hv.accuracy, 
      hv.dims = hv.dims
    )
  
  })

  results[parallel$batch] <- res

  catn("Finishing up.")

  stopCluster(parallel$cl)
  
  prev_highest_it <- as.integer(readLines(parallel$highest.iteration))
  
  if (length(prev_highest_it) == 0) prev_highest_it <- 0

  highest_iteration <- max(unlist(lapply(results, function(res) res$iteration)))
  
  if (highest_iteration >= prev_highest_it) {
    catn("Highest iteration:", highcat(highest_iteration))
    catn("Expected total iteration:", highcat(length(spec.list)))
    writeLines(as.character(highest_iteration), parallel$highest.iteration)
  }
  
  end_timer(parallell_timer)

  vebcat("Hypervolume sequence completed succesfully", color = "seqSuccess")
}
