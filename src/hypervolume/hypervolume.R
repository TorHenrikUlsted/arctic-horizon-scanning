source_all("./src/hypervolume/components")

process_species <- function(spec.dt, spec.name, method, points.projection = "longlat", verbose = F, iteration, warn.file, err.file) {
  vebcat("Initiating data processing protocol.", color = "funInit")

  biovars_world <- rast(".outputs/setup/region/biovars-world-subset.tif")

  vebcat("Starting species data process.", veb = verbose)

  withCallingHandlers(
    {
      sp_points <- prepare_species(
        spec.dt = spec.dt,
        projection = "longlat",
        verbose = verbose
      )
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when preparing species", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when preparing species", iteration = iteration)
  )

  if (is.null(sp_points)) {
    catn("Excluding", spec.name, "from further processing.")
    return(list(
      excluded = TRUE
    ))
  }

  vebprint(sp_points, verbose, "Processed environment data sample:")
  vebcat("Using biovars:", highcat(names(biovars_world)), veb = verbose)
  vebcat("Starting processing of environment data.", veb = verbose)

  withCallingHandlers(
    {
      sp_mat <- prepare_environment(
        sp_points,
        biovars_world,
        verbose = verbose
      )
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when preparing environment", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when preparing environment", iteration = iteration)
  )

  vebprint(head(sp_mat, 3), verbose, "Processed environment data sample:")

  vebcat("Data processing protocol completed successfully.", color = "funSuccess")

  return(sp_mat)
}

hv_analysis <- function(spec.mat, method, spec.name, incl_threshold, accuracy, iteration, hv.dir, cores.max.high, verbose, warn.file, err.file) {
  vebcat("Initiating hypervolume sequence", color = "funInit")

  locks_dir <- paste0(hv.dir, "/locks/hypervolume")
  proj_dir <- paste0(hv.dir, "/projection/", spec.name)

  create_dir_if(c(locks_dir, proj_dir))

  catn("Reading raster files.")
  biovars_region <- rast(".outputs/setup/region/biovars-region-subset.tif")

  # check if need to setup, or read setup region
  region_hv <- setup_hv_region(
    biovars_region,
    out.dir = ".outputs/setup/region",
    method = method
  )

  # Remove those with too few observations
  # Determined by the log(observations) being less than number of dimensions
  catn("Analyzing for curse of dimensionality. \n")

  withCallingHandlers(
    {
      nobs <- nrow(spec.mat)
      ndim <- ncol(spec.mat)
      spp <- ceiling((10^(3 + sqrt(ndim))) / nobs)
      nrp <- spp * nobs

      cat(sprintf("%15s | %15s | %20s \n", "n_observations", "n_dimensions", "log(n_observations)"))
      cat(sprintf("%15.4f | %15.4f | %20.4f \n", nobs, ndim, log(nobs)))
      catn("Samples per point:", highcat(spp))

      if (log(nobs) <= ndim) {
        catn("log(n_observations)", highcat(log(nobs)), "smaller than n_dimensions,", highcat(ndim), "excluding from further analysis. \n")

        return(list(
          n_observations = nobs,
          n_dimensions = ndim,
          samples_per_point = spp,
          random_points = nrp,
          excluded = TRUE,
          analyzed_hv_stats = c(rep(0, 2), rep(0, 2)),
          included_sp = c(rep(TRUE, 0), rep(FALSE, nobs))
        ))
      }
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when analyzing curse of dimensionality", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when analyzing curse of dimensionality", iteration = iteration)
  )

  ## Inclusion test to eliminate obvious non-overlaps
  catn("Computing inclusion analysis. \n")

  withCallingHandlers(
    {
      included_sp <- hypervolume_inclusion_test(
        region_hv,
        sp_mat,
        reduction.factor = 1,
        fast.or.accurate = accuracy,
        fast.method.distance.factor = 1,
        accurate.method.threshold = quantile(region_hv@ValueAtRandomPoints, 0.5),
        verbose = verbose
      )

      vebcat("Inclusion analysis completed successfully", color = "funSuccess")

      catn(
        "Number of TRUE / FALSE = OVERLAP values:",
        highcat(sum(included_sp == T)),
        "/",
        highcat(sum(included_sp == F)),
        "=",
        highcat(format(sum(included_sp == T) / length(included_sp), nsmall = 2, big.mark = ","))
      )

      if (any(included_sp == T)) {
        vebcat("Included for further hypervolume analysis.", color = "proSuccess")
      } else {
        vebcat("Excluded from further hypervolume analysis.", color = "nonFatalError")

        return(list(
          n_observations = nobs,
          n_dimensions = ndim,
          samples_per_point = spp,
          random_points = nrp,
          excluded = TRUE,
          analyzed_hv_stats = c(rep(0, 2), rep(0, 2)),
          included_sp = included_sp
        ))
      }
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when analyzing hypervolume inclusion", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when analyzing hypervolume inclusion", iteration = iteration)
  )

  analysis_msg <- FALSE
  mm <- FALSE
  ml <- get_mem_usage("total", "gb") * 0.8
  # Wait for a spot in the queue system
  while (TRUE) {
    mc <- get_mem_usage("used", "gb")

    if (!analysis_msg) {
      catn("The node has entered the hypervolume analysis queue...")
      analysis_msg <- TRUE
    }

    if (mc >= ml) {
      if (!mm) {
        catn("The node waiting for more available RAM before initiating the hypervolume analysis.")
        mm <- TRUE
      }

      Sys.sleep(10)
    } else {
      if (is.locked(locks_dir, lock.n = cores.max.high)) {
        Sys.sleep(10)
      } else {
        Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
        lock_analysis <- lock(locks_dir, lock.n = cores.max.high)
        if (!is.null(lock_analysis) && file.exists(lock_analysis)) {
          break
        }
      }
    }
  }

  analysis_msg <- FALSE

  catn("The node has exited the queue.")

  ## If included, continue with hypervolume analysis
  catn("Computing hypervolume analysis.")
  withCallingHandlers(
    {
      sp_hv <- hypervolume(spec.mat, name = spec.name, method = method, verbose = T)
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when analyzing hypervolume", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when analyzing hypervolume", iteration = iteration)
  )

  withCallingHandlers(
    {
      analyzed_hv_stats <- analyze_hv_stats(region_hv, sp_hv, spec.name, verbose = verbose)

      vebprint(analyzed_hv_stats, text = "Hypervolume Statistics:")

      # If the region is unique and does not overlap with any of the species
      if (analyze_hv_stats[[4]] == 1) {
        catn(sp_hv@Name, colcat("Excluded from further hypervolume analysis.", color = "nonFatalError"))

        vebcat("Unlocking analysis lock.", veb = verbose)
        unlock(lock_analysis)

        return(list(
          n_observations = nobs,
          n_dimensions = ndim,
          samples_per_point = spp,
          random_points = nrp,
          excluded = TRUE,
          analyzed_hv_stats = analyzed_hv_stats,
          included_sp = included_sp
        ))
      } else {
        catn(sp_hv@Name, colcat("Included for further hypervolume analysis.", color = "proSuccess"))
      }
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when analyzing statistics", iteration = j),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when analyzing statistics", iteration = j)
  )

  withCallingHandlers(
    {
      catn("Projecting inclusion analysis.")
      catn("Using the", accuracy, "inclusion accuracy.")
      ## Projections for inclusion
      for (threshold in incl_threshold) {
        inc_project <- hypervolume_project(
          sp_hv,
          biovars_region,
          type = "inclusion",
          fast.or.accurate = accuracy,
          accurate.method.threshold = quantile(sp_hv@ValueAtRandomPoints, threshold),
          verbose = TRUE
        )

        names(inc_project) <- ("inclusionScore")
        out_file <- paste0(proj_dir, "/inclusion-", threshold, ".tif")

        vebcat("Writing out raster file:", colcat(out_file), color = "output")
        writeRaster(inc_project, out_file, overwrite = TRUE)

        rm(inc_project)
        invisible(gc())
      }

      catn("Projecting probability analysis.")
      ## Projections for probability
      prob_proj <- hypervolume_project(
        sp_hv,
        biovars_region,
        type = "probability",
        verbose = verbose
      )

      names(prob_proj) <- ("suitabilityScore")
      out_file <- paste0(proj_dir, "/probability.tif")

      vebcat("Writing out raster file:", colcat(out_file), color = "output")
      writeRaster(prob_proj, out_file, overwrite = T)
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when projecting hypervolume", iteration = iteration),
    error = function(e) err(e, err.file = err.file, err.txt = "Error when projecting hypervolume", iteration = iteration)
  )

  vebcat("Unlocking analysis lock.", veb = verbose)
  unlock(lock_analysis)

  rm(prob_proj)
  invisible(gc())

  veb("Hypervolume sequence completed successfully.", color = "funSuccess")

  return(list(
    n_observations = nobs,
    n_dimensions = ndim,
    samples_per_point = spp,
    random_points = nrp,
    excluded = FALSE,
    analyzed_hv_stats = analyzed_hv_stats,
    included_sp = included_sp
  ))
}
