source_all("./R/hypervolume/components")

process_species <- function(spec.dt, spec.name, process.dir, method, points.projection = "longlat", verbose = F, iteration, warn.file, err.file) {
  # Create a list to track objects that need cleanup
  objects_to_clean <- list()

  tryCatch(
    {
      vebcat("Initiating data processing protocol.", color = "funInit")

      biovars_world <- rast(paste0(build_climate_path(), "/biovars-world-subset.tif"))

      vebcat("Starting species data process.", veb = verbose)

      withCallingHandlers(
        {
          sp_points <- prepare_species(
            dt = spec.dt,
            process.dir = process.dir,
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

      vebprint(sp_points, verbose, "Processed species data sample:")
      vebcat("Using biovars:", names(biovars_world), veb = verbose)
      vebcat("Starting processing of environment data.", veb = verbose)

      withCallingHandlers(
        {
          sp_mat <- prepare_environment(
            sp_points,
            biovars_world,
            verbose = verbose
          )

          # remove sp_points to save ram
          rm(sp_points)
          invisible(gc())
        },
        warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when preparing environment", iteration = iteration),
        error = function(e) err(e, err.file = err.file, err.txt = "Error when preparing environment", iteration = iteration)
      )

      vebprint(head(sp_mat, 3), verbose, "Processed environment data sample:")

      vebcat("Data processing protocol completed successfully.", color = "funSuccess")

      return(sp_mat)
    },
    error = function(e) {
      # Handle any errors that weren't caught by withCallingHandlers
      err(e, err.file = err.file, err.txt = "Unexpected error in process_species", iteration = iteration)
    },
    finally = function() {
      # Cleanup code that runs regardless of success or failure
      tryCatch(
        {
          # Clean up tracked objects
          for (obj in objects_to_clean) {
            if (!is.null(obj)) {
              if (inherits(obj, "SpatRaster") || inherits(obj, "SpatVector")) {
                terra::tmpFiles(remove = TRUE) # Clean terra temporary files
              }
              rm(obj)
            }
          }
          # Clean up environment and force garbage collection
          rm(list = ls(environment()))
          gc(full = TRUE, reset = TRUE)
        },
        error = function(e) {
          warning("Error during cleanup of process_species: ", e$message)
        }
      )
    }
  )
}

hv_analysis <- function(spec.mat, method, spec.name, incl_threshold, accuracy, iteration, hv.dir, lock.dir, proj.dir, cores.max.high, verbose, warn.file, err.file) {
  sp_hv <- NULL
  region_hv <- NULL
  lock_analysis <- NULL
  start_mem <- get_mem_usage("used", "gb")

  tryCatch(
    {
      vebcat("Initiating hypervolume sequence", color = "funInit")
      proj_dir <- paste0(proj.dir, "/", spec.name)
      create_dir_if(lock.dir)

      catn("Reading raster files.")
      biovars_region <- rast(paste0(build_climate_path(), "/biovars-region-subset.tif"))

      # check if need to setup, or read setup region
      region_hv <- setup_hv_region(
        biovars_region,
        out.dir = build_climate_path(),
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
              analyzed_hv_stats = c(0, 0, 1, 1),
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
            spec.mat,
            reduction.factor = 1,
            fast.or.accurate = accuracy,
            fast.method.distance.factor = 1,
            accurate.method.threshold = quantile(region_hv@ValueAtRandomPoints, 0.5),
            verbose = verbose
          )

          vebcat("Inclusion analysis completed successfully", color = "funSuccess")

          catn(
            "Number of TRUE / FALSE = OVERLAP values:",
            highcat(sum(included_sp == TRUE)),
            "/",
            highcat(sum(included_sp == FALSE)),
            "=",
            highcat(format(sum(included_sp == TRUE) / length(included_sp), nsmall = 2, big.mark = ","))
          )

          if (any(included_sp == TRUE)) {
            vebcat("Included for further hypervolume analysis.", color = "proSuccess")
          } else {
            vebcat("Excluded from further hypervolume analysis.", color = "nonFatalError")

            return(list(
              n_observations = nobs,
              n_dimensions = ndim,
              samples_per_point = spp,
              random_points = nrp,
              excluded = TRUE,
              analyzed_hv_stats = c(0, 0, 1, 1),
              included_sp = included_sp
            ))
          }
        },
        warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when analyzing hypervolume inclusion", iteration = iteration),
        error = function(e) err(e, err.file = err.file, err.txt = "Error when analyzing hypervolume inclusion", iteration = iteration)
      )

      # Wait for memory and a spot in the queue system
      repeat {
        # Check memory first
        if (!mem_check(
          identifier = paste0("analysis_", spec.name),
          ram.use = paste0(hv.dir, "/logs/ram-usage.txt"),
          custom_msg = paste("Waiting in hypervolume analysis queue for", spec.name),
          verbose = verbose
        )) {
          # If memory is OK, check for lock availability
          if (!is.locked(lock.dir, lock.n = cores.max.high)) {
            # Try to acquire lock
            Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
            lock_analysis <- lock(lock.dir, lock.n = cores.max.high, paste0("Locked by ", iteration, "_", spec.name))
            if (!is.null(lock_analysis) && file.exists(lock_analysis)) {
              break # Exit loop when memory is OK and lock is acquired
            }
          }
          Sys.sleep(10) # Wait before trying again if lock not acquired
        }
      }

      catn("The node has exited the queue.")

      ## If included, continue with hypervolume analysis
      catn("Computing hypervolume analysis.")
      withCallingHandlers(
        {
          sp_hv <- hypervolume(spec.mat, name = spec.name, method = method, verbose = T)
        },
        warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when computing hypervolume", iteration = iteration),
        error = function(e) {
          unlock(lock_analysis)
          err(e, err.file = err.file, err.txt = "Error when computing hypervolume", iteration = iteration)
        }
      )

      catn("Analyzing statistics.")
      withCallingHandlers(
        {
          # Add check for valid hypervolume objects
          if (is.null(sp_hv) || is.null(region_hv)) {
            stop("Invalid hypervolume objects - one or both are NULL")
          }

          analyzed_hv_stats <- analyze_hv_stats(sp_hv, region_hv, spec.name, verbose = verbose)

          vebprint(analyzed_hv_stats, text = "Hypervolume Statistics:")

          # If the region is unique and does not overlap with any of the species
          if (analyzed_hv_stats[[4]] == 1) {
            catn(spec.name, colcat("Excluded from further hypervolume analysis.", color = "nonFatalError"))

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
            catn(spec.name, colcat("Included for further hypervolume analysis.", color = "proSuccess"))
          }
        },
        warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when analyzing statistics", iteration = iteration),
        error = function(e) {
          unlock(lock_analysis)
          err(e, err.file = err.file, err.txt = "Error when analyzing statistics", iteration = iteration)
        }
      )

      withCallingHandlers(
        {
          create_dir_if(proj_dir)
          catn("Projecting inclusion analysis.")
          catn("Using the", paste0("'", accuracy, "'"), "inclusion accuracy.")
          catn("Using", incl_threshold, "threshold(s).")
          ## Projections for inclusion
          for (threshold in incl_threshold) {
            out_file <- paste0(proj_dir, "/", threshold, "-inclusion", ".tif")
            if (file.exists(out_file)) {
              catn(threshold, "-inclusion analysis file already exists.")
              next
            }

            inc_project <- hypervolume_project(
              sp_hv,
              biovars_region,
              type = "inclusion",
              fast.or.accurate = accuracy,
              accurate.method.threshold = quantile(sp_hv@ValueAtRandomPoints, threshold),
              verbose = verbose
            )

            names(inc_project) <- ("inclusionScore")

            inc_project <- check_crs(
              inc_project,
              projection = config$simulation$projection,
              projection.method = "near",
              res = config$projection$raster_scale_m,
              verbose = verbose
            )

            vebcat("Writing out raster file:", colcat(out_file), color = "output")
            writeRaster(inc_project, out_file, overwrite = TRUE)

            rm(inc_project)
            invisible(gc())
          }


          out_file <- paste0(proj_dir, "/probability.tif")

          if (!file.exists(out_file)) {
            catn("Projecting probability analysis.")

            prob_project <- hypervolume_project(
              sp_hv,
              biovars_region,
              type = "probability",
              verbose = verbose
            )

            names(prob_project) <- ("suitabilityScore")

            prob_project <- check_crs(
              prob_project,
              projection = config$simulation$projection,
              projection.method = "bilinear",
              res = config$projection$raster_scale_m,
              verbose = verbose
            )

            vebcat("Writing out raster file:", colcat(out_file), color = "output")
            writeRaster(prob_project, out_file, overwrite = TRUE)

            rm(prob_project)
            invisible(gc())
          } else {
            catn("Probability analysis file already exists.")
          }
        },
        warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when projecting hypervolume", iteration = iteration),
        error = function(e) {
          unlock(lock_analysis)
          err(e, err.file = err.file, err.txt = "Error when projecting hypervolume", iteration = iteration)
        }
      )

      rm(biovars_region)
      invisible(gc())

      vebcat("Unlocking analysis lock.", veb = verbose)
      unlock(lock_analysis)

      vebcat("Hypervolume analysis completed successfully.", color = "funSuccess")

      return(list(
        n_observations = nobs,
        n_dimensions = ndim,
        samples_per_point = spp,
        random_points = nrp,
        excluded = FALSE,
        analyzed_hv_stats = analyzed_hv_stats,
        included_sp = included_sp
      ))
    },
    error = function(e) {
      # If there's an error, make sure to unlock before propagating
      if (!is.null(lock_analysis) && file.exists(lock_analysis)) {
        unlock(lock_analysis)
      }
      stop(e)
    },
    finally = function() {
      # Clean up large objects first
      if (exists("sp_hv")) {
        rm(sp_hv)
        gc(full = TRUE)
      }
      if (exists("region_hv")) {
        rm(region_hv)
        gc(full = TRUE)
      }
      if (exists("biovars_region")) {
        rm(biovars_region)
        gc(full = TRUE)
      }

      # Clean up any remaining lock
      if (!is.null(lock_analysis) && file.exists(lock_analysis)) {
        unlock(lock_analysis)
      }

      # Clean up everything else
      rm(list = ls(environment()))
      gc(full = TRUE)

      # Log memory usage
      end_mem <- get_mem_usage("used", "gb")
      mem_diff <- end_mem - start_mem
      if (mem_diff > 1) { # If we retained more than 1GB
        warning("Memory not fully cleaned up in hv_analysis for ", spec.name, "in iteration", iteration, ". Diff: ", mem_diff, "GB")
      }
    }
  )
}
