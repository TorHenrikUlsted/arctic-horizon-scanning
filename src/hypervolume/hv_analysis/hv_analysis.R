hv_analysis <- function(sp_mat, biovars_region, region_hv, method, spec.name, proj.incl.t, accuracy, hv.projection, verbose, iteration, proj.dir, lock_hv_dir, cores.max.high, warn, err) {
  cat(blue("Initiating hypervolume sequence \n"))
  
  ## Inclusion test to eliminate obvious non-overlaps
  cat("Computing inclusion analysis. \n")
  
  withCallingHandlers(
    {
      nobs <- nrow(sp_mat)
      ndim <- ncol(sp_mat)
      spp <- ceiling((10^(3 + sqrt(ndim))) / nobs)
      nrp <- spp * nobs
      
      cat(sprintf("%15s | %15s | %20s \n", "n_observations", "n_dimensions", "log(n_observations)"))
      cat(sprintf("%15.4f | %15.4f | %20.4f \n", nobs, ndim, log(nobs)))
      cat("Samples per point:", cc$lightSteelBlue(spp), "\n")
      
      if (log(nobs) <= ndim) {
        cat("log(n_observations)", cc$lightSteelBlue(log(nobs)), "smaller than n_dimensions,", cc$lightSteelBlue(ndim), "excluding from further analysis. \n")
        
        return(list(
          n_observations = nobs,
          n_dimensions = ndim,
          samples_per_point = spp,
          random_points = nrp,
          excluded = T,
          analyzed_hv_stats = c(rep(0, 2), rep(0, 2)),
          included_sp = c(rep(TRUE, 0), rep(FALSE, nobs))
        ))
      }
      
      included_sp <- hypervolume_inclusion_test(
        region_hv,
        sp_mat,
        reduction.factor = 1,
        fast.or.accurate = "accurate",
        fast.method.distance.factor = 1,
        accurate.method.threshold = quantile(region_hv@ValueAtRandomPoints, 0.5),
        verbose = verbose
      )
      
      invisible(gc())
      
      cat(cc$lightGreen("inclusion analysis completed successfully \n"))
      
      cat("Number of TRUE / FALSE = OVERLAP values:", cc$lightSteelBlue(sum(included_sp == T)), "/", cc$lightSteelBlue(sum(included_sp == F)), "=", cc$lightSteelBlue(format(sum(included_sp == T) / length(included_sp), nsmall = 2, big.mark = ",")), "\n")
      if (any(included_sp == T)) {
        cat(green("Included for further hypervolume analysis. \n"))
      } else {
        cat(red("Excluded from further hypervolume analysis. \n"))
        
        return(list(
          n_observations = nobs,
          n_dimensions = ndim,
          samples_per_point = spp,
          random_points = nrp,
          excluded = T,
          analyzed_hv_stats = c(rep(0, 2), rep(0, 2)),
          included_sp = included_sp
        ))
      }
    },
    warning = function(w) warn(w, warn_txt = "Warning when analyzing hypervolume inclusion in iteration"),
    error = function(e) err(e, err_txt = "Error when analyzing hypervolume inclusion in iteration")
  )
  
  invisible(gc())
  
  analysis_msg <- FALSE
  mm <- FALSE
  ml <- get_mem_usage("total", "gb") * 0.8
  # Wait for a spot in the queue system
  while (TRUE) {
    mc <- get_mem_usage("used", "gb")
    
    if (!analysis_msg) {
      cat("The node has entered the hypervolume analysis queue... \n")
      analysis_msg <- TRUE
    }
    
    if (mc >= ml) {
      if (!mm) {
        cat("The node waiting for more available RAM before initiating the hypervolume analysis. \n")
        mm <- TRUE
      }
      
      Sys.sleep(10)
      
    } else {
      
      if (is.locked(lock_hv_dir, lock.n = cores.max.high)) {
        Sys.sleep(10)
      } else {
        Sys.sleep(runif(1, 0, 1))  # Add a random delay between 0 and 1 second
        lock_analysis <- lock(lock_hv_dir, lock.n = cores.max.high)
        if (!is.null(lock_analysis) && file.exists(lock_analysis)) {
          break
        }
      }
    }
  }
  
  analysis_msg <- FALSE
  
  cat("The node has exited the queue. \n")
  
  ## If included, continue with hypervolume analysis
  cat("Computing hypervolume analysis. \n")
  withCallingHandlers(
    {
      sp_hv <- hypervolume(sp_mat, name = spec.name, method = method, verbose = T)
      
    },
    warning = function(w) warn(w, warn_txt = "Warning when analyzing hypervolume in iteration"),
    error = function(e) err(e, err_txt = "Error when analyzing hypervolume in iteration")
  )
  
  invisible(gc())
  
  withCallingHandlers(
    {
      analyzed_hv_stats <- analyze_hv_stats(region_hv, sp_hv, spec.name, verbose = T)
      
      cat("Hypervolume Statistics", method, "\n")
      print(analyzed_hv_stats)
      
      if (1 - analyzed_hv_stats[[4]] > 0) {
        cat(sp_hv@Name, green("Included for further hypervolume analysis. \n"))
      } else {
        cat(sp_hv@Name, red("Excluded from further hypervolume analysis. \n"))
        
        if (verbose) cat("Unlocking analysis lock. \n")
        unlock(lock_analysis)
        
        return(list(
          n_observations = nobs,
          n_dimensions = ndim,
          samples_per_point = spp,
          random_points = nrp,
          excluded = T,
          analyzed_hv_stats = analyzed_hv_stats,
          included_sp = included_sp
        ))
      }
    },
    warning = function(w) warn(w, warn_txt = "Warning when analyzing statistics in iteration"),
    error = function(e) err(e, err_txt = "Error when analyzing statistics in iteration")
  )
  
  invisible(gc())
  
  withCallingHandlers(
    {
      proj_dir <- paste0(proj.dir, "/", method, "/", spec.name)
      create_dir_if(proj_dir)
      
      cat("Projecting inclusion analysis. \n")
      cat("Using the", accuracy, "inclusion accuracy. \n")
      ## Projections for inclusion
      for (threshold in proj.incl.t) {
        inc_project <- hypervolume_project(
          sp_hv,
          biovars_region,
          type = "inclusion",
          fast.or.accurate = accuracy,
          accurate.method.threshold = quantile(sp_hv@ValueAtRandomPoints, threshold),
          verbose = T
        )
        names(inc_project) <- ("inclusionScore")
        if (hv.projection == "laea") {
          cat("Reprojecting to laea. \n")
          inc_proj <- terra::project(inc_project, crs(laea_crs), method="near")
        } else cat("Keeping longlat projection. \n")
        writeRaster(inc_project, paste0(proj_dir, "/inclusion-", threshold, ".tif"), overwrite = T)
        rm(inc_project)
        invisible(gc())
      }
      
      
      cat("Projecting probability analysis. \n")
      ## Projections for probability
      prob_proj <- hypervolume_project(sp_hv, biovars_region, type = "probability", verbose = T)
      
      names(prob_proj) <- ("suitabilityScore")
      
      if (hv.projection == "laea") {
        cat("Reprojecting to laea. \n")
        prob_proj <- terra::project(prob_proj, crs(laea_crs), method = "bilinear")
      } else cat("Keeping the longlat projection. \n")
      
      if (verbose) cat("Writing out raster file:", paste0(proj_dir, "/probability.tif"), "\n")
      writeRaster(prob_proj, paste0(proj_dir, "/probability.tif"), overwrite = T)
    },
    warning = function(w) warn(w, warn_txt = "Warning when projecting hypervolume in iteration"),
    error = function(e) err(e, err_txt = "Error when projecting hypervolume in iteration")
  )
  
  if (verbose) cat("Unlocking analysis lock. \n")
  unlock(lock_analysis)
  
  rm(prob_proj)
  invisible(gc())
  
  cat(cc$lightGreen("Hypervolume sequence completed successfully. \n\n"))
  
  return(list(
    n_observations = nobs,
    n_dimensions = ndim,
    samples_per_point = spp,
    random_points = nrp,
    excluded = F,
    analyzed_hv_stats = analyzed_hv_stats,
    included_sp = included_sp
  ))
}