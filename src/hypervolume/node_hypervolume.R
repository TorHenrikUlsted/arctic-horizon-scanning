node_hypervolume <- function(
    process.dir,
    iteration,
    spec.list,
    columns.to.read,
    min.disk.space,
    cores.max.high,
    verbose = FALSE,
    hv.incl.threshold, 
    hv.method, 
    hv.accuracy, 
    hv.dims
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
    fun = function(spec, spec.name) {
      
      skip_to_end <- FALSE
      proj_dir <- paste0(process.dir, "/projections")
      pro_locks_dir <- paste0(node$locks, "/hv-locks")
      
      create_dir_if(c(proj_dir, pro_locks_dir))
      
      catn("Log observations:", log(nobs))
      catn("Expected dimensions:", length(hv.dims))
      
      if (log(nobs) < length(hv.dims)) {
        catn("Number of observations are less than the expected dimensions.")
        skip_to_end <- TRUE
      }
      
      while (!skip_to_end) {
        catn("Running main sequence.")
        
        processed_data <- process_species(
          spec.dt = spec, 
          spec.name = spec.name, 
          method = hv.method, 
          iteration = iteration,
          verbose = verbose,
          warn.file = node$warn,
          err.file = node$err
        )
        
        if (is.list(processed_data) && processed_data$excluded == TRUE) {
          skip_to_end <- TRUE
          break
        }
        
        analyzed_hv <- hv_analysis(
          spec.mat = processed_data,
          method = hv.method, 
          spec.name = spec.name, 
          incl_threshold = hv.incl.threshold, 
          accuracy = hv.accuracy, 
          projection = hv.projection, 
          iteration = iteration, 
          hv.dir = process.dir,
          lock.dir = lock_hv_dir, 
          cores.max.high = cores.max.high,
          verbose = verbose,
          warn.file = node$warn,
          err.file = node$err
        )
        
        catn("Appending data to csv file.")
        
        final_res <- data.table(
          scientificName = spec.name,
          iteration = iteration,
          observations = analyzed_hv[[1]],
          dimensions = analyzed_hv[[2]],
          samplesPerPoint = analyzed_hv[[3]],
          randomPoints = analyzed_hv[[4]],
          excluded = analyzed_hv[[5]],
          jaccard = analyzed_hv[[6]][[1]],
          sørensen = analyzed_hv[[6]][[2]],
          fracVolumeSpecies = analyzed_hv[[6]][[3]],
          fracVolumeRegion = analyzed_hv[[6]][[4]],
          overlapRegion = 1 - analyzed_hv[[6]][[4]],
          includedOverlap = sum(analyzed_hv[[7]] == T) / length(analyzed_hv[[7]])
        )
        # If we want the number of true and false values we can do: round(observations * includedOverlap) = true values, and: observations - true values = false values
        
        break
      }
      
        if (skip_to_end) {
          catn("Skipping to end.")
          catn("Appending data to csv file.")
          
          final_res <- data.table(
            scientificName = spec.name,
            iteration = iteration,
            observations = nobs,
            dimensions = length(hv.dims),
            samplesPerPoint = 0,
            randomPoints = 0,
            excluded = TRUE,
            jaccard = 0,
            sørensen = 0,
            fracVolumeSpecies = 0,
            fracVolumeRegion = 0,
            overlapRegion = 1,
            includedOverlap = 0
          )
        }
        
        fwrite(final_res, paste0(process.dir, "/stats.csv"), append = T, bom = T)
        
        vebprint(final_res, text = "Appended data:")
      }
    )

  rm(list = setdiff(ls(), "iteration"))

  invisible(gc())

  list(
    iteration = iteration
  )
}
