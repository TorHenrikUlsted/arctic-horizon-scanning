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
    hv.incl.threshold) {
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
    level3Name = character(0),
    level3Code = character(0),
    level2Code = integer(0),
    level1Code = integer(0),
    level3Long = numeric(0),
    level3Lat = numeric(0)
  )

  cols_to_select <- c(
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
  )

  if (!file.exists(stats_file)) {
    fwrite(init_dt, stats_file, row.names = F, bom = T)
  }

  custom_exports <- c(
    hv_dir,
    hv.method,
    hv.accuracy,
    hv.dims,
    hv.incl.threshold,
    cols_to_select
  )

  custom_evals <- c(
    "./R/setup/setup_sequence.R",
    "./R/hypervolume/hypervolume.R",
    "./R/hypervolume/node_hypervolume.R",
    "./R/hypervolume/parallel_hypervolume.R"
  )

  if (!file.exists(removed_file)) {
    spec_count_dt <- exclude_observations(
      spec.list = spec.list,
      dimensions = hv.dims,
      dt.construct = init_dt,
      cols.select = cols_to_select,
      out.file = removed_file,
      method = "median",
      suppress = TRUE,
      verbose = verbose
    )

    # spec_removed <- spec_count_dt[removed == TRUE, ]
    # spec_removed[, filename := NULL]
    #
    # catn(highcat(nrow(spec_removed)), "Species with too few observations were removed.")
    #
    # spec_count_dt <- spec_count_dt[removed == FALSE, ]
  }

  if (!file.exists(spec_list_file)) {
    spec_list <- unlist(optimize_queue(
      spec_count_dt,
      cores.max,
      high_ram_threshold = 0.2,
      verbose = FALSE
    ))

    writeLines(unlist(spec_list), spec_list_file)
    catn("Adding to markdown file")
    mdwrite(
      config$files$post_seq_md,
      text = paste0(
        "1;Hypervolume Sequence\n\n",
        "Species removed before analysis because of too few occurrences: ",
        "**", nrow(spec_removed), "**  ",
        "\n\nSpecies input into the hypervolume sequence: ",
        "**", length(spec_list), "**"
      ),
    )
  } else {
    spec_list <- readLines(spec_list_file)
    catn("Read from file with", highcat(length(spec_list)), "species")
  }

  if (is.character(iterations)) {
    spec_list <- unlist(spec_list)
    spec_list <- spec_list[grep(iterations, gsub(config$species$file_separator, " ", spec_list))]
    iterations <- length(spec_list)
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
    return(vebcat("Hypervolume already finished a run for this list.", color = "funSuccess"))
  }

  vebprint(clusterEvalQ(parallel$cl, ls())[[1]], veb = verbose, text = "cluster variables sample:")

  # ETC calculation
  hv_setup_time <- paste0("./outputs/setup/hypervolume/", hv.method, "-sequence/setup-check/avg-time.txt")
  if (file.exists(hv_setup_time)) {
    base_time_iteration <- as.numeric(readLines(hv_setup_time)) * 60
  } else {
    base_time_iteration <- 25 * 60
  }

  etc <- calculate_etc(
    timer.res = base_time_iteration,
    cores = parallel$cores,
    data.length = length(parallel$batch)
  )

  tryCatch(
    {
      clusterApplyLB(parallel$cl, parallel$batch, function(j) {
        initial_objects <- ls(envir = .GlobalEnv, all.names = TRUE)

        tryCatch({
          repeat {
            if (!mem_check(
              identifier = paste("node", j),
              ram.use = parallel$ram.use,
              verbose = verbose
            )) {
              break # Exit loop when memory is OK
            }
          }

          node_hypervolume(
            process.dir = par.dir, # par.dir because it comes from inside the parllel setup
            iteration = j,
            spec.list = spec_list,
            columns.to.read = cols_to_select,
            min.disk.space = min.disk.space,
            cores.max.high = cores.max.high,
            init.dt = init_dt,
            verbose = verbose,
            hv.incl.threshold = hv.incl.threshold,
            hv.method = hv.method,
            hv.accuracy = hv.accuracy,
            hv.dims = hv.dims
          )
        }, error = function(e) {
          err_msg <- paste("Error in node_hypervolume in iteration", j, ":", e$message)
          warning(err_msg)
        }, finally = {
          # Clean up worker's own environment
          rm(list = setdiff(ls(all.names = TRUE), initial_objects))
          gc(full = TRUE)
        })
      })
    },
    error = function(e) {
      vebcat("Error in hypervolume parallel process ~ stopping clusters and closing all connections.", color = "fatalError")
      stop(e)
    }, finally = {
      catn("Cleaning up parallel process.")
      cleanup_files(paste0(hv_dir, "/locks"), recursive = TRUE, verbose = verbose)
      stopCluster(parallel$cl)
      closeAllConnections()
    }
  )

  catn("Finishing up")

  highest_iteration <- as.integer(readLines(parallel$highest.iteration))
  do_not_return <- FALSE

  if (highest_iteration >= length(spec_list)) {
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

  rm(list = ls(environment()))
  gc(full = TRUE)
}
