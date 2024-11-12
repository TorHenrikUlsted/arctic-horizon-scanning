node_hypervolume <- function(
    process.dir,
    iteration,
    spec.list,
    columns.to.read,
    min.disk.space = 1,
    cores.max.high = 1,
    init.dt,
    verbose = FALSE,
    hv.incl.threshold,
    hv.method,
    hv.accuracy,
    hv.dims) {
  
  tryCatch({
    node <- setup_node(
      pro.dir = process.dir,
      iteration = iteration,
      min.disk.space = min.disk.space,
      verbose = verbose
    )
  }, error = function(e) {
    vebcat("Error when setting up node in iteration:", iteration, color = "fatalError")
    stop(e)
  })
  
  
  tryCatch(
    {
      process_node(
        pro.dir = process.dir,
        identifier = node$identifier,
        lock.dir = node$locks,
        lock.setup = node$lock.setup,
        node.it.log = node$it.log,
        node.init.log = node$init.log,
        iteration = iteration,
        spec.list = spec.list,
        columns.to.read = columns.to.read,
        init.dt = init.dt,
        hv.incl.threshold = hv.incl.threshold,
        verbose = verbose,
        warn = node$warn.file,
        err = node$err.file,
        fun = function(spec, spec.name) {
          on.exit({
            # Clean up any resources created within this function
            rm(list = ls(environment()))
            gc(full = TRUE)
          })
          
          skip_to_end <- FALSE
          proj_dir <- paste0(process.dir, "/projections")
          pro_locks_dir <- paste0(dirname(node$locks), "/hypervolume")

          create_dir_if(proj_dir, pro_locks_dir)

          nobs <- nrow(spec)

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
              process.dir = process.dir,
              method = hv.method,
              iteration = iteration,
              verbose = verbose,
              warn.file = node$warn,
              err.file = node$err
            )

            if (is.list(processed_data) && processed_data$excluded == TRUE) {
              vebcat("Setting skip_to_end to TRUE.", veb = verbose)
              skip_to_end <- TRUE
              break
            }

            analyzed_hv <- hv_analysis(
              spec.mat = processed_data,
              method = hv.method,
              spec.name = spec.name,
              incl_threshold = hv.incl.threshold,
              accuracy = hv.accuracy,
              iteration = iteration,
              hv.dir = process.dir,
              lock.dir = pro_locks_dir,
              proj.dir = proj_dir,
              cores.max.high = cores.max.high,
              verbose = verbose,
              warn.file = node$warn,
              err.file = node$err
            )

            catn("Appending data to csv file.")

            final_res <- data.table(
              cleanName = gsub(config$species$file_separator, " ", spec.name),
              iteration = iteration,
              observations = analyzed_hv[[1]],
              dimensions = analyzed_hv[[2]],
              samplesPerPoint = analyzed_hv[[3]],
              randomPoints = analyzed_hv[[4]],
              excluded = analyzed_hv[[5]],
              jaccard = analyzed_hv[[6]][[1]],
              sorensen = analyzed_hv[[6]][[2]],
              fracVolumeSpecies = analyzed_hv[[6]][[3]],
              fracVolumeRegion = analyzed_hv[[6]][[4]],
              realizedNiche = 1 - analyzed_hv[[6]][[3]],
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
              cleanName = gsub(config$species$file_separator, " ", spec.name),
              iteration = iteration,
              observations = nobs,
              dimensions = length(hv.dims),
              samplesPerPoint = 0,
              randomPoints = 0,
              excluded = TRUE,
              jaccard = 0,
              sorensen = 0,
              fracVolumeSpecies = 1,
              fracVolumeRegion = 1,
              realizedNiche = 0,
              overlapRegion = 0,
              includedOverlap = 0
            )
          }

          vebprint(final_res, text = "Hypervolume data return:")

          return(final_res)
        }
      )

    },
    warning = function(w) warn(w, warn.file = node$warn.file, warn.txt = "Warning in node_hypervolume", iteration = iteration),
    error = function(e) {
      err(e, err.file = node$err.file, err.txt = "Error in node_hypervolume", iteration = iteration)
      stop(e)
    }, finally = function() {
      catn("Cleaning up node environment")
      closeAllConnections()
      invisible(gc(full = TRUE))
    }
  )
  
  
  
  return(catn("Node returning."))
}
