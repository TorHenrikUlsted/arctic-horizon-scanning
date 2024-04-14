hypervolume_sequence <- function(
    spec.list,
    iterations = NULL, 
    cores.max.high = 1,
    cores.max = 1,
    min.disk.space = 1,
    verbose = FALSE,
    hv.method, 
    hv.accuracy, 
    hv.dims, 
    hv.incl.threshold
    ) {
  on.exit(closeAllConnections())
  
  vebcat("Initiating hypervolume sequence", color = "seqInit")
  
  parallell_timer <- start_timer("parallell_timer")
  
  process_dir <- paste0("./outputs/hypervolume/", method, "-sequence")
  
  stats_file <- paste0(hv_dir, "/stats.csv")
  
  if (!file.exists(stats_file)) {
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
      mean_long = numeric(0),
      mean_lat = numeric(0)
    )
    
    fwrite(init_dt, stats_file, row.names = F, bom = T)
  }
  
  custom_exports = c(
    hv.incl.threshold, 
    hv.method, 
    hv.accuracy, 
    hv.dims
  )
  
  custom_evals = c(
    "./src/utils/utils.R",
    "./src/setup/setup.R",
    "./src/hypervolume/data_acquisition/data_acquisition.R",
    "./src/hypervolume/data_analysis/data_analysis.R",
    "./src/hypervolume/data_processing/data_processing.R",
    "./src/hypervolume/hv_analysis/hv_analysis.R",
    "./src/hypervolume/node.R"
  )
  
  parallel <- setup_parallel(
    par.dir = process_dir,
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

  res <- clusterApplyLB(cl, parallel$batch, function(j) {
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
    
    node_processing(
      process.dir = hv.dir,
      iteration = j,
      min.disk.space = min.disk.space,
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
      hv.incl.threshold = incl.threshold, 
      hv.method = , 
      hv.accuracy, 
      hv.dims
    )
  
  })

  results[parallel$batch] <- res

  catn("Finishing up.")

  stopCluster(cl)
  
  prev_highest_it <- as.integer(readLines(parallel$highest.iteration))

  highest_iteration <- max(unlist(lapply(results, function(res) res$iteration)))
  
  if (highest_iteration > prev_highest_it) writeLines(as.character(highest_iteration), highest_it_file)
  
  end_timer(parallell_timer)

  veb("Hypervolume sequence completed succesfully", color = "seqSuccess")
}
