hv_analysis <- function(sp_mat, biovars_region, region_hv, method, spec.name, proj.incl.t, accuracy, hv.projection, verbose, iteration, proj.dir, lock_hv_dir, cores.max.high, warn, err) {
  vebcat("Initiating hypervolume sequence", color = "funInit")
  
  ## Inclusion test to eliminate obvious non-overlaps
  catn("Computing inclusion analysis. \n")
  
  withCallingHandlers(
    {
      nobs <- nrow(sp_mat)
      ndim <- ncol(sp_mat)
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
      
      included_sp <- hypervolume_inclusion_test(
        region_hv,
        sp_mat,
        reduction.factor = 1,
        fast.or.accurate = "accurate",
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
    warning = function(w) warn(w, warn_txt = "Warning when analyzing hypervolume inclusion in iteration"),
    error = function(e) err(e, err_txt = "Error when analyzing hypervolume inclusion in iteration")
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
  
  catn("The node has exited the queue.")
  
  ## If included, continue with hypervolume analysis
  catn("Computing hypervolume analysis.")
  withCallingHandlers(
    {
      sp_hv <- hypervolume(sp_mat, name = spec.name, method = method, verbose = T)
      
    },
    warning = function(w) warn(w, warn_txt = "Warning when analyzing hypervolume in iteration"),
    error = function(e) err(e, err_txt = "Error when analyzing hypervolume in iteration")
  )
  
  withCallingHandlers(
    {
      analyzed_hv_stats <- analyze_hv_stats(region_hv, sp_hv, spec.name, verbose = T)
      
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
    warning = function(w) warn(w, warn_txt = "Warning when analyzing statistics in iteration"),
    error = function(e) err(e, err_txt = "Error when analyzing statistics in iteration")
  )
  
  withCallingHandlers(
    {
      proj_dir <- paste0(proj.dir, "/", method, "/", spec.name)
      create_dir_if(proj_dir)
      
      catn("Projecting inclusion analysis.")
      catn("Using the", accuracy, "inclusion accuracy.")
      ## Projections for inclusion
      for (threshold in proj.incl.t) {
        inc_project <- hypervolume_project(
          sp_hv,
          biovars_region,
          type = "inclusion",
          fast.or.accurate = accuracy,
          accurate.method.threshold = quantile(sp_hv@ValueAtRandomPoints, threshold),
          verbose = TRUE
        )
        
        names(inc_project) <- ("inclusionScore")
        
        if (hv.projection == "laea") {
          catn("Reprojecting to laea.")
          inc_proj <- terra::project(inc_project, crs(laea_crs), method="near")
        } else catn("Keeping longlat projection.")
        
        out_file <- paste0(proj_dir, "/inclusion-", threshold, ".tif")
        vebcat("Writing out raster file:", colcat(out_file), color = "output")
        writeRaster(inc_project, out_file, overwrite = TRUE)
        
        rm(inc_project)
        invisible(gc())
      }
      
      
      catn("Projecting probability analysis.")
      ## Projections for probability
      prob_proj <- hypervolume_project(sp_hv, biovars_region, type = "probability", verbose = T)
      
      names(prob_proj) <- ("suitabilityScore")
      
      if (hv.projection == "laea") {
        catn("Reprojecting to laea.")
        prob_proj <- terra::project(prob_proj, crs(laea_crs), method = "bilinear")
      } else catn("Keeping the longlat projection.")
      
      out_file <- paste0(proj_dir, "/probability.tif")
      vebcat("Writing out raster file:", colcat(out_file), color = "output")
      writeRaster(prob_proj, out_file, overwrite = T)
    },
    warning = function(w) warn(w, warn_txt = "Warning when projecting hypervolume in iteration"),
    error = function(e) err(e, err_txt = "Error when projecting hypervolume in iteration")
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