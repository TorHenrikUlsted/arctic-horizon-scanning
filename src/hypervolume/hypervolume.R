source("./src/hypervolume/data_acquisition/acquire_data.R")
source("./src/hypervolume/data_analysis/analyze_data.R")
source("./src/hypervolume/data_processing/process_data.R")

###########################################
#                                         #
# ---------- Data Acquisition ----------- #
#                                         #
###########################################

data_acquisition <- function(show.plot, method, verbose, iteration, warn, err) {
  cat(blue("Initiating data acquisition protocol \n"))

  if (verbose) cat(blue("Acquiring regions. \n"))

  shapefiles <- c(
    cavm = "./resources/region/cavm-noice/cavm-noice.shp"
  )
  withCallingHandlers(
    {
      regions <- import_regions(shapefiles, "./outputs/hypervolume/data_acquisition/region")
    }, 
    warning = function(w) warn(w, warn_txt = "Warning when importing regions in iteration"),
    error = function(e) err(e, err_txt = "Error when importing regions in iteration")
  )
  
  cavm_floreg <- terra::split(regions$cavm, regions$cavm$FLOREG)

  for (i in seq_along(cavm_floreg)) {
    names(cavm_floreg)[i] <- paste("floreg", unique(cavm_floreg[[i]]$FLOREG), sep = "_")
    if (verbose) cat("Renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(names(cavm_floreg)[i]), "\n")
  }

  if (verbose) cat(cc$lightGreen("Region acquisition completed successfully \n"))

  if (verbose) cat(blue("Acquiring WorldClim biovariables. \n"))

  withCallingHandlers(
    {
      biovars_world <- get_wc_data(show.plot = show.plot, verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when getting worldClim data in iteration"),
    error = function(e) err(e, err_txt = "Error when getting worldClim data in iteration")
  )

  if (verbose) cat(cc$lightGreen("Biovars_world acquired successfully \n"))

  if (verbose) cat(blue("Scaling biovars_world \n"))

  withCallingHandlers(
    {
      biovars_world <- scale_biovars(biovars_world, verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when scaling biovars_world in iteration"),
    error = function(e) err(e, err_txt = "Error when scaling biovars_world in iteration")
  )

  if (verbose) cat(cc$lightGreen("Biovars_world scaled successfully \n"))

  if (verbose) cat(blue("Acquiring biovars_region \n"))

  withCallingHandlers(
    {
      biovars_region <- acquire_region_data(biovars_world, regions, projection = "longlat", verbose = verbose)
    },
    warning =  function(w) warn(w, warn_txt = "Warning when acquiring region data in iteration"),
    error = function(e) err(e, err_txt = "Error when acquiring region data in iteration")
  )

  if (verbose) cat(cc$lightGreen("Biovars_region acquired successfully \n"))

  if (verbose) cat(cc$lightGreen("Acquiring biovars for floristic regions \n"))

  withCallingHandlers(
    {
      biovars_floreg <- acquire_region_data(biovars_world, cavm_floreg, projection = "longlat", verbose = verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when acquiring floristic region data in iteration"),
    error = function(e) err(e, err_txt = "Error when when acquiring floristic region data in iteration")
  )

  if (verbose) cat(cc$lightGreen("Biovars_floreg acquired successfully \n"))

  cat(cc$lightGreen("Data acquisition protocol completed successfully \n"))

  return(list(
    biovars_world,
    biovars_region,
    biovars_floreg
  ))
}

###########################################
#                                         #
# ------------ Data Analysis ------------ #
#                                         #
###########################################

data_analysis <- function(biovars_world, biovars_region, biovars_floreg, method, show.plot, verbose, iteration, warn, err) {
  cat(blue("Initiating data analysis protocol \n"))
  withCallingHandlers(
    {
      analyzed_data <- analyze_correlation(biovars_region, file.out = "./outputs/hypervolume/data_analysis/correlation", threshold = 0.5)
    },
    warning = function(w) warn(w, warn_txt = "Warning when analyzing correlation in iteration"),
    error = function(e) err(e, err_txt = "Error when analyzing correlation in iteration")
  )

  withCallingHandlers(
    {
      # choose wanted correlation dimensions
      biovars_world <- terra::subset(biovars_world, c(18, 10, 3, 4))
      biovars_region <- terra::subset(biovars_region, c(18, 10, 3, 4))
      for (i in seq_along(biovars_floreg)) {
        subset <- biovars_floreg[[i]][[c(18, 10, 3, 4)]]

        biovars_floreg[[i]] <- subset
      }
    },
    warning = function(w) warn(w, warn_txt = "Warning when subsetting regions in iteration"),
    error = function(e) err(e, err_txt = "Error when when subsetting regions in iteration")
  )

  withCallingHandlers(
    {
      region_hv <- setup_region_hv(biovars_region, "cavm", method = method)
    },
    warning = function(w) warn(w, warn_txt = "Warning when setting up region hypervolume in iteration"),
    error = function(e) err(e, err_txt = "Error when setting up region hypervolume in iteration")
  )

  cat(cc$lightGreen("Data analysis protocol completed successfully \n"))

  return(list(
    biovars_world,
    biovars_region,
    biovars_floreg,
    region_hv
  ))
}

###########################################
#                                         #
# ----------- Data Processing ----------- #
#                                         #
###########################################

data_processing <- function(sp_df, biovars_world, spec.name, method, projection = "longlat", verbose = F, iteration, warn, err) {
  cat(blue("Initiating data processing protocol \n"))

  if (verbose) cat(blue("Processing species data.\n"))

  withCallingHandlers(
    {
      sp_points <- prepare_species(sp_df, projection = "longlat", verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when preparing species in iteration"),
    error = function(e) err(e, err_txt = "Error when preparing species in iteration")
  )

  if (verbose) {
    cat("Processed environment data sample: \n")
    print(sp_points)
    cat(cc$lightGreen("Species preperation completed successfully. \n"))
  }

  if (verbose) {
    cat(blue("Processing environment data.\n"))
    cat("Using biovars:", cc$lightSteelBlue(names(biovars_world)), "\n")
  }

  withCallingHandlers(
    {
      sp_mat <- prepare_environment(sp_points, biovars_world, verbose)
    },
    warning = function(w) warn(w, warn_txt = "Warning when preparing environment in iteration"),
    error = function(e) err(e, err_txt = "Error when preparing environment in iteration")
  )

  if (verbose) {
    cat(cc$lightGreen("Environment preperation completed successfully. \n"))
    cat("Processed environment data sample: \n")
    print(head(sp_mat, 3))
  }


  cat(cc$lightGreen("Data processing protocol completed successfully. \n"))

  return(sp_mat)
}

###########################################
#                                         #
# ------------- Hypervolume ------------- #
#                                         #
###########################################


hv_analysis <- function(sp_mat, biovars_region, region_hv, method, spec.name, proj.incl.t, accuracy, project, verbose, iteration, warn, err) {
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
      proj_dir <- paste0("./outputs/hypervolume/projections/", method, "/", gsub(" ", "-", spec.name))
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
         if (project == "laea") {
           cat("Reprojecting to laea. \n")
           inc_proj <- terra::project(inc_project, crs(laea_crs))
         } else cat("Keeping longlat projection. \n")
         writeRaster(inc_project, paste0(proj_dir, "/inclusion-", threshold, ".tif"), overwrite = T)
         rm(inc_project)
         invisible(gc())
       }


       cat("Projecting probability analysis. \n")
       ## Projections for probability
       prob_proj <- hypervolume_project(sp_hv, biovars_region, type = "probability", verbose = T)
      
       names(prob_proj) <- ("suitabilityScore")
      
       if (project == "laea") {
         cat("Reprojecting to laea. \n")
         prob_proj <- terra::project(prob_proj, crs(laea_crs))
       } else cat("Keeping the longlat projection. \n")
      
       writeRaster(prob_proj, paste0(proj_dir, "/probability.tif"), overwrite = T)
    },
    warning = function(w) warn(w, warn_txt = "Warning when projecting hypervolume in iteration"),
    error = function(e) err(e, err_txt = "Error when projecting hypervolume in iteration")
  )

  rm(prob_proj)
  invisible(gc())
  
  cat(cc$lightGreen("Hypervolume sequence completed successfully. \n"))

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
