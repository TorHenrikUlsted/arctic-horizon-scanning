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
  stats_dir <- paste0(hv_dir, "/stats")
  
  create_dir_if(stats_dir)
  
  stats_file <- paste0(stats_dir, "/stats.csv")
  removed_file <- paste0(stats_dir, "/removed-species.csv")
  spec_list_file <- paste0(stats_dir, "/spec-iteration-list.txt")
  
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
    potentialRealizedNiche = numeric(0),
    potentialOverlapRegion = numeric(0),
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
  
  if (!file.exists(spec_list_file)) {
    spec_count_dt <- count_observations(
      spec.list = spec.list,
      dimensions = hv.dims,
      method = "median",
      verbose = verbose
    )
    
    spec_removed <- spec_count_dt[removed == TRUE, ]
    spec_removed[, filename := NULL]
    
    catn(highcat(nrow(spec_removed)), "Species with too few observations were removed.")
    
    fwrite(spec_removed, removed_file, bom = TRUE)
    
    spec_count_dt <- spec_count_dt[removed == FALSE, ]
    
    setorder(spec_count_dt, medianLat)
    
    spec_list <- spec_count_dt$filename
    
    vebprint(head(spec_list, 10), verbose, "species list sample:")
    
    writeLines(spec_list, spec_list_file)
    
    mdwrite(
      post_seq_nums,
      heading = paste0("1;Hypervolume Sequence\n\n",
                       "Species removed before analysis because of too few occurrences: ",
                       "**",nrow(spec_removed),"**  ",
                       "Species input into the hypervolume sequence: ",
                       "**",length(spec_list),"**"
                      ),
    )
  } else {
    spec_list <- readLines(spec_list_file)
  }
  
  catn("Starting hypervolume with", highcat(length(spec_list)), "species.")
  
  parallel <- setup_parallel(
    par.dir = hv_dir,
    spec.list = spec_list,
    iterations = iterations,
    cores.max = cores.max,
    cores.max.high = cores.max.high,
    min.disk.space = min.disk.space,
    verbose = verbose,
    custom.exports = custom_exports,
    custom.evals = custom_evals
  )
  
  if (parallel$finished) {
    return(catn("Hypervolume already finished a run for this list."))
  }
  
  vebprint(clusterEvalQ(parallel$cl, ls()), veb = verbose, text = "All cluster variables:")

 tryCatch({
   res <- clusterApplyLB(parallel$cl, parallel$batch, function(j) {
     invisible(gc())
     
     ram_msg <- FALSE
     # RAM check
     mem_used_gb <- get_mem_usage(type = "used", format = "gb")
     mem_limit_gb <- mem_limit / 1024^3
     
     while (mem_used_gb >= mem_limit_gb & !file.exists(paste0(par.dir, "/logs/escape.txt"))) {
       if (!ram_msg) {
         ram_con <- file(parallel$ram.use, open = "a")
         writeLines(paste0("RAM usage ", mem_used_gb, " is above the maximum ", mem_limit_gb, " Waiting with node", j), ram_con)
         close(ram_con)
         ram_msg = TRUE
       }
       Sys.sleep(60)  # Wait for 60 seconds before checking again
       Sys.sleep(runif(1, 0, 1)) # Add random seconds between 0 and 1 to apply difference if multiple nodes are in queue
       mem_used_gb <- get_mem_usage(type = "used", format = "gb")
     }
     
     node_hypervolume(
       process.dir = par.dir, # par.dir because it comes from inside the parllel setup
       iteration = j,
       spec.list = spec_list,
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
 },
  error = function(e) {
    stopCluster(parallel$cl)
    closeAllConnections()
    vebcat("Error in hypervolume parallel process, stopping clusters and closing all connections.", color = "fatalError")
    print(e)
  }
)

  catn("Finishing up.")

  stopCluster(parallel$cl)
  
  highest_iteration <- as.integer(readLines(parallel$highest.iteration))
  do_not_return <- FALSE
  
  if (highest_iteration == length(spec_list)) {
    vebcat("Parallel process finished all iterations.", color = "funSuccess")
  } else {
    do_not_return <- TRUE
    vebcat("Parallel process did not finish all iterations.", color = "nonFatalError")
  }
  catn("Highest iteration/Expected highest iteration:")
  vebcat(highest_iteration, "/", length(spec_list), color = "indicator")
  
  # Check if all nodes ran
  stats_file_its <- fread(paste0(stats_dir, "/stats.csv"), select = "iteration")
  stats_file_its <- unique(stats_file_its$iteration)
  expected_its <- 1:(length(spec_list))
  
  missing_iterations <- setdiff(expected_its, stats_file_its)
    
  if (length(missing_iterations) > 0) {
    do_not_return <- TRUE
    vebcat("Missing iterations found:\n", missing_iterations, color = "fatalError")
  } 
  
  end_timer(parallell_timer)

  vebcat("Hypervolume sequence completed succesfully", color = "seqSuccess")
  
  if (do_not_return) stop("Stopping process from going to visualization sequence because of missing data.")
}
