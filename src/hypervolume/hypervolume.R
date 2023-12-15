# source("./src/utils/utils.R")
#source("./src/setup/setup.R")
source("./src/hypervolume/data_acquisition/acquire_data.R")
source("./src/hypervolume/data_analysis/analyze_data.R")
source("./src/hypervolume/data_processing/process_data.R")

###########################################
#                                         #
# ---------- Data Acquisition ----------- #
#                                         #
###########################################

data_acquisition <- function(show_plot, verbose,  iteration) {
  cat(blue("Initiating data acquisition protocol \n"))

  if (verbose) cat(blue("Acquiring regions. \n"))

  shapefiles <- c(
    cavm = "./resources/region/cavm_noice/cavm_noice.shp"
  )
  withCallingHandlers(
    {
      regions <- import_regions(shapefiles, "./outputs/data_acquisition/region/logs")
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when importing regions in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when importing regions in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
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
      biovars_world <- get_wc_data(show_plot = show_plot, verbose = F)
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when acquiring WorldClim data in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when acquiring WorldClim data in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
  )
  if (verbose) cat(cc$lightGreen("Biovars_world acquired successfully \n"))

  if (verbose) cat(blue("Scaling biovars_world \n"))
  withCallingHandlers(
    {
      biovars_world <- scale_biovars(biovars_world, verbose = verbose)
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when scaling biovars_world in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when caling biovars_world in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
  )
  if (verbose) cat(cc$lightGreen("Biovars_world scaled successfully \n"))

  if (verbose) cat(blue("Acquiring biovars_region \n"))
  withCallingHandlers(
    {
      biovars_region <- acquire_region_data(biovars_world, regions, projection = "longlat", verbose = verbose)
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when acquiring region data in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when region data in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
  )
  if (verbose) cat(cc$lightGreen("Biovars_region acquired successfully \n"))

  if (verbose) cat(cc$lightGreen("Acquiring biovars for floristic regions \n"))
  withCallingHandlers(
    {
      biovars_floreg <- acquire_region_data(biovars_world, cavm_floreg, projection = "longlat", verbose = verbose)
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when acquiring floristic region data in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when acquiring floristic region data in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
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

data_analysis <- function(biovars_world, biovars_region, biovars_floreg, show_plot, verbose, iteration) {
  cat(blue("Initiating data analysis protocol \n"))
  withCallingHandlers(
    {
      # analyzed_data <- analyze_correlation(biovars_region, threshold = 0.5)
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when analysing correlation data in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when analysing correlation data in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
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
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when subsetting biovars in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when subsetting biovars in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
  )
  withCallingHandlers(
    {
      if (file.exists("./outputs/data_analysis/hypervolume/region/cavm/hypervolume_box.rds")) {
        cat("Region hypervolume already exists. \n")
        region_hv <- readRDS("./outputs/data_analysis/hypervolume/region/cavm/hypervolume_box.rds")
      } else {
        region_hv_log <- "./outputs/data_analysis/hypervolume/region/box_messages.txt"

        try(region_hv_log <- file(region_hv_log, open = "at"))

        sink(region_hv_log, type = "message")

        region_hv <- analyze_region_hv(biovars_region, "cavm", method = "box", samples.per.point = 1, verbose = T)

        sink(type = "message")

        close(region_hv_log)
      }
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when analysing hypervolume for region in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when analysing hypervolume for region in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
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

data_processing <- function(sp_df, biovars_world, spec.name, projection = "longlat", verbose = F, iteration) {
  cat(blue("Initiating data processing protocol \n"))

  if (verbose) cat(blue("Processing species data.\n"))
  withCallingHandlers(
    {
      sp_points <- prepare_species(sp_df, projection, verbose)
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when preparing species in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when preparing species in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
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
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning when preparing environment in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error when preparing environment in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
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


hv_analysis <- function(sp_mat, biovars_region, region_hv, method, spec.name, verbose, iteration) {
  cat(blue("Initiating hypervolume sequence \n"))

  ## Inclusion test to eliminate obvious non-overlaps
  cat("Computing inclusion analysis. \n")

  withCallingHandlers(
    {
      nobs <- nrow(sp_mat)
      ndim <- ncol(sp_mat)
      spp <- ceiling((10^(3 + sqrt(ndim)))/nobs)
      nrp <- spp * nobs
      
      cat("Number of observations:", cc$lightSteelBlue(nobs), "\n")
      cat("Samples per point:", cc$lightSteelBlue(spp), "\n")
      
      included_sp <- hypervolume_inclusion_test(region_hv,
        sp_mat,
        reduction.factor = 1,
        fast.or.accurate = "accurate",
        fast.method.distance.factor = 1,
        accurate.method.threshold = quantile(region_hv@ValueAtRandomPoints, 0.5),
        verbose = verbose
      )
      cat(cc$lightGreen("inclusion analysis completed successfully \n"))

      cat("Number of TRUE / FALSE = OVERLAP values:", cc$lightSteelBlue(sum(included_sp == T)), "/", cc$lightSteelBlue(sum(included_sp == F)), "=", cc$lightSteelBlue(format(sum(included_sp == T) / length(included_sp), nsmall = 2, big.mark = ",")), "\n")
      if (any(included_sp == T)) {
        if (nobs <= ndim) {
          cat("Number of observations", cc$lightSteelBlue(nobs), "smaller than number of dimensions", cc$lightSteelBlue(ndim), "excluding from further analysis. \n")
          return(list(
            included_sp = included_sp, 
            analyzed_hv_stats = rep(0, 4),
            random_points = nrp,
            samples_per_point = spp,
            n_observations = nobs,
            n_dimensions = ndim
          ))
        }
        
        cat(green("Included for further hypervolume analysis. \n"))
        
      } else {
        cat(red("Excluded from further hypervolume analysis. \n"))
        return(list(
          included_sp = included_sp, 
          analyzed_hv_stats = rep(0, 4),
          random_points = nrp,
          samples_per_point = spp,
          n_observations = nobs,
          n_dimensions = ndim
        ))
      }
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning in inclusion analysis in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error in inclusion analysis in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
  )

  ## If included, continue with hypervolume analysis
  cat("Computing hypervolume analysis. \n")
  withCallingHandlers(
    {
      sp_hv <- hypervolume(sp_mat, name = spec.name, method = method, verbose = T)
      
    },
    warning = function(w) {
      warn_msg <- conditionMessage(w)
      sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
      cat("Warning ins pecies hypervolume analysis in iteration", iteration, ":", warn_msg, "\n")
      sink()
      invokeRestart(findRestart('muffleWarning'))
    },
    error = function(e) {
      err_msg <- conditionMessage(e)
      sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
      cat("Error in species hypervolume analysis in iteration", iteration, ":", err_msg, "\n")
      sink()
      stop(e)
    }
  )
  
  withCallingHandlers(
  {
    analyzed_hv_stats <- analyze_hv_stats(region_hv, sp_hv, spec.name, verbose = T)

    print(analyzed_hv_stats)

    if (1 - analyzed_hv_stats[[4]] > 0) {
      cat(sp_hv@Name, green("Included for further hypervolume analysis. \n"))
    } else {
      cat(sp_hv@Name, red("Excluded from further hypervolume analysis. \n"))

      return(list(
        included_sp = included_sp,
        analyzed_hv_stats = analyzed_hv_stats,
        random_points = nrp,
        samples_per_point = spp,
        n_observations = nobs,
        n_dimensions = ndim
      ))
    }
  },
  warning = function(w) {
    warn_msg <- conditionMessage(w)
    sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
    cat("Warning in species hypervolume statistics analysis in iteration", iteration, ":", warn_msg, "\n")
    sink()
    invokeRestart(findRestart('muffleWarning'))
  },
  error = function(e) {
    err_msg <- conditionMessage(e)
    sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
    cat("Error in species hypervolume statistics analysis in iteration", iteration, ":", err_msg, "\n")
    sink()
    stop(e)
  }
)
    
  # withCallingHandlers(
  #   {
  #     # Move onto projections
  #     proj_dir <- paste0("./outputs/hypervolume/projections/", method, "/", gsub(" ", "-", spec.name))
  #     create_dir_if(proj_dir)
  # 
  #     cat("Projecting inclusion analysis. \n")
  #     ## Projections for inclusion
  #     sp_inc_project <- hypervolume_project(sp_hv, biovars_region, type = "inclusion", verbose = T)
  # 
  #     names(sp_inc_project) <- ("suitabilityScore")
  # 
  #     laea_inc_proj <- terra::project(sp_inc_project, crs(laea_crs))
  # 
  #     writeRaster(laea_inc_proj, paste0(proj_dir, "/inclusion.tif"), overwrite = T)
  # 
  #     cat("Projecting probability analysis. \n")
  #     ## Projections for probability
  #     sp_prob_project <- hypervolume_project(sp_hv, biovars_region, type = "probability", verbose = T)
  # 
  #     names(sp_prob_project) <- ("suitabilityScore")
  # 
  #     laea_prob_proj <- terra::project(sp_prob_project, crs(laea_crs))
  # 
  #     writeRaster(laea_prob_proj, paste0(proj_dir, "/probability.tif"), overwrite = T)
  #   },
  #   warning = function(w) {
  #     warn_msg <- conditionMessage(w)
  #     sink(paste0("./outputs/hypervolume/logs/", method, "-warning.txt"), append = T)
  #     cat("Warning in inclusion analysis in iteration", iteration, ":", warn_msg, "\n")
  #     sink()
  #     invokeRestart(findRestart('muffleWarning'))
  #   },
  #   error = function(e) {
  #     err_msg <- conditionMessage(e)
  #     sink(paste0("./outputs/hypervolume/logs/", method, "-error.txt"), append = T)
  #     cat("Error in inclusion analysis in iteration", iteration, ":", err_msg, "\n")
  #     sink()
  #     stop(e)
  #   }
  # )


  cat(cc$lightGreen("Hypervolume sequence completed successfully. \n"))
  return(list(
    included_sp = included_sp,
    analyzed_hv_stats = analyzed_hv_stats,
    random_points = nrp,
    samples_per_point = spp,
    n_observations = nobs,
    n_dimensions = ndim
  ))
}
