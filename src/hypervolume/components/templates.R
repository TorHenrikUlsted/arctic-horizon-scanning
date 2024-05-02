parallel_sequence <- function(
    spec.list,
    iterations = NULL, 
    cores.max.high = 1,
    cores.max = 1,
    min.disk.space,
    cluster.evals,
    verbose = FALSE#,
    # ---- Custom parameters ---- #
) {
  on.exit(closeAllConnections())
  
  vebcat("Initiating hypervolume sequence", color = "seqInit")
  
  parallell_timer <- start_timer("parallell_timer")
  
# ---- Custom Block ---- #
  
  # Initialize things here
  process_dir <- paste0("./outputs/somedir")
  
  custom_exports = c(
    ""
  )
  
  custom_evals = c(
    ""
  )
  
# ---- Custom Block ---- #
  
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
        "species", 
        "decimalLongitude", 
        "decimalLatitude", 
        "coordinateUncertaintyInMeters", 
        "coordinatePrecision", 
        "countryCode", 
        "stateProvince", 
        "year"
      ),
      verbose = verbose#,
      # ---- Custom parameters ---- #
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

#################
# Node template
#################

# This template will create directories based off of process.dir, and log files will initiate
node_process <- function(
    process.dir,
    iteration,
    min.disk.space,
    spec.list,
    columns.to.read,
    verbose = FALSE#,
    # Your own parameters
) {
  
  node <- setup_node(
    pro.dir = process.dir,
    iteration = iteration,
    min.disk.space = min.disk.space,
    verbose = verbose
  )
  
  process_node(
    pro.dir = process.dir,
    iteration = iteration,
    identifier = node$identifier,
    spec.list = spec.list,
    columns.to.read = columns.to.read,
    verbose = verbose,
    fun = function() {
      
    }
  )
  
  rm(list = setdiff(ls(), "j"))
  
  invisible(gc())
  
  list(
    iteration = j
  )
}