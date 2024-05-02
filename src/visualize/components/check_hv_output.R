check_hv_output <- function(spec.list, hv.dir, hv.method, hv.projection, hv.inc.t = 0.5, verbose = F) {
  vebcat("Checking Hypervolume output", color = "funInit")
  
  proj_dir <- paste0(hv.dir, "/projections")
  stats_dir <- paste0(hv.dir, "/stats")
  logs_dir <- "./outputs/visualize/logs"

  spec_list <- gsub("\\.csv$", "", spec.list)
  spec_list_dt <- data.table(species = basename(spec_list), dirName = spec_list)
  spec_list_dt$iteration <- seq_len(nrow(spec_list_dt))

  proj_dirs <- list.dirs(path = proj_dir, recursive = FALSE, full.names = TRUE)

  catn("Running projection check.")
  for (i in seq_along(proj_dirs)) {
    method_dir <- proj_dirs[[i]]

    missing_rast_log <- paste0(logs_dir, "/", basename(method_dir), "-missing-rasters.csv")
    create_file_if(missing_rast_log)

    wrong_crs_log <- paste0(logs_dir, "/", basename(method_dir), "-wrong_crs.csv")
    create_file_if(wrong_crs_log)

    wrong_ext_log <- paste0(logs_dir, "/", basename(method_dir), "-wrong_exent.csv")
    create_file_if(wrong_ext_log)

    sp_dirs <- list.dirs(path = method_dir, full.names = TRUE)
    sp_dirs_dt <- data.table(species = basename(sp_dirs), dirName = sp_dirs, iteration = NA)
    sp_dirs_dt <- sp_dirs_dt[-1, ]

    if (any(duplicated(sp_dirs_dt$species))) {
      dup_sp <- duplicated(sp_dirs_dt$species)
      dup_names <- sp_dirs_dt$species[sup_sp]
      catn("Found", catcol(sum(dup_sp), "highlight"), "duplicated species directories")
      print(dup_names)
      stop("Needs manual check")
    } else {
      vebcat("Could not find any projection directory duplicates.", color = "proSuccess")
    }

    merged_dt <- merge(sp_dirs_dt, spec_list_dt, by = "species")
    merged_dt <- merged_dt[, -c("dirName.y", "iteration.x")]
    names(merged_dt) <- c("species", "dirName", "hvIteration")

    missing_rast <- data.table(species = character(), dirName = character(), hvIteration = numeric())
    wrong_crs <- data.table(species = character(), hvIteration = numeric(), crs = character(), fileName = character())
    wrong_extent <- data.table(species = character(), hvIteration = numeric(), ext = character(), expectedExt = character(), res = numeric(), expectedRes = numeric(), fileName = character())
    catn("Expected crs:", crs(hv.projection, proj = TRUE))
    # Inside the species directory of the Hypervolume output
    for (j in 1:nrow(merged_dt)) {
      tryCatch(
        {
          sp <- merged_dt$species[j]
          sp_dir <- merged_dt$dirName[j]
          hv_it_num <- merged_dt$hvIteration[j]

          vebcat("Species:", sp, veb = verbose)
          vebcat("Directory:", sp_dir, veb = verbose)
          vebcat("Iteration number:", hv_it_num, veb = verbose)

          rast_files <- list.files(path = sp_dir, pattern = "\\.tif$", full.names = TRUE)

          # Check number of files within the directory
          cat("\rChecking directory", j, "/", nrow(sp_dirs_dt))
          flush.console()

          ## If the dir has less files than expected, add the iteration of the species to a file

          if (length(rast_files) < (length(hv.inc.t) + 1)) {
            missing_rast_i <- data.table(species = sp, dirName = sp_dir, hvIteration = hv_it_num)
            missing_rast <- rbind(missing_rast, missing_rast_i)
          }

          # Check each file in the directory
          for (k in seq_along(rast_files)) {
            rast_file <- rast_files[[k]]

            if (!file.exists(rast_file)) {
              vebcat("\nCannot find file:", rast_file, color = "nonFatalError", veb = verbose)
            } else {
              vebcat("\nFound file:", rast_file, color = "proSuccess", veb = verbose)
            }

            rast <- terra::rast(rast_file)

            vebcat("Rast:", veb = verbose)
            vebprint(rast, veb = verbose)

            rast_ext <- terra::ext(rast)

            if (basename(rast_file) == "inclusion-0.5.tif") {
              # Compare extents
              if (!all(ext(-180, 180, 55.7916666666667, 83.625) == rast_ext)) {
                vebcat("Raster has mismatched extent:", color = "nonFatalError", veb = verbose)
                vebcat("basename:", basename(rast_file), veb = verbose)
                vebcat("Extent:", as.character(rast_ext), "/", ext(-180, 180, 55.7916666666667, 83.625), veb = verbose)

                wrong_ext_i <- data.table(
                  species = sp,
                  hvIteration = hv_it_num,
                  ext = as.character(rast_ext),
                  expectedExt = "-180, 180, 55.7916666666667, 83.625",
                  res = res(rast),
                  expectedRes = 0.04166667,
                  fileName = rast_file
                )
                wrong_extent <- rbind(wrong_extent, wrong_ext_i)
              }
            }

            if (basename(rast_file) == "probability.tif") {
              # Compare extents
              if (!all(ext(-180, 180, 40.795614586063, 90) == rast_ext)) {
                vebcat("Raster has mismatched extent:", color = "nonFatalError", veb = verbose)
                vebcat("basename:", basename(rast_file), veb = verbose)
                vebcat("Extent:", as.character(rast_ext), "/", ext(-180, 180, 40.795614586063, 90), veb = verbose)

                wrong_ext_i <- data.table(
                  species = sp,
                  hvIteration = hv_it_num,
                  ext = as.character(rast_ext),
                  expectedExt = "-180, 180, 40.795614586063, 90",
                  res = res(rast),
                  expectedRes = c(0.02076963, 0.02077011),
                  fileName = rast_file
                )
                wrong_extent <- rbind(wrong_extent, wrong_ext_i)
              }
            }

            raster_crs <- terra::crs(rast, proj = TRUE)
            
            vebcat("raster_crs:", veb = verbose)
            vebprint(raster_crs, veb = verbose)

            if (identical(raster_crs, crs(hv.projection, proj = TRUE))) {
              vebcat("Raster has correct projection.", color = "proSuccess", veb = verbose)
              vebcat("CRS:", raster_crs, "/", crs(hv.projection, proj = TRUE), veb = verbose)
            } else if (!identical(raster_crs, crs(hv.projection, proj = TRUE))) {
              if (verbose) {
                vebcat("Raster has incorrect projection.", color = "nonFatalError", veb = verbose)
                vebcat("CRS:", raster_crs, "/", crs(hv.projection, proj = TRUE), veb = verbose)
              }

              wrong_crs_i <- data.table(species = sp, hvIteration = hv_it_num, crs = as.character(raster_crs), fileName = rast_file)
              wrong_crs <- rbind(wrong_crs, wrong_crs_i)
            }
          }
        },
        error = function(e) {
          print(paste("Error at iteration", j, ":", e$message))
        }
      )
    }; catn()

    if (nrow(missing_rast) > 0 || nrow(wrong_crs) > 0 || nrow(wrong_extent) > 0) {
      if (nrow(missing_rast) > 0) {
        vebcat("hypervolume output failed the projection check with missing rasters.", color = "nonFatalError")
        
        fwrite(missing_rast, missing_rast_log, bom = TRUE)
        
        catn("Missing rasters:", highcat(nrow(missing_rast)), "/", highcat(nrow(sp_dirs_dt) * (length(hv.inc.t) + 1)))
      }

      if (nrow(wrong_crs) > 0) {
        vebcat("hypervolume output failed the projection check with wrong crs.", color = "nonFatalError")
        
        fwrite(wrong_crs, wrong_crs_log, bom = TRUE)
        
        catn("Rasters with wrong crs:", highcat(nrow(wrong_crs)), "/", highcat(nrow(sp_dirs_dt) * (length(hv.inc.t) + 1)))
      }

      if (nrow(wrong_extent) > 0) {
        vebcat("hypervolume output failed the projection check with wrong extents.", color = "nonFatalError")
        
        fwrite(wrong_extent, wrong_ext_log, bom = TRUE)
        
        catn("Rasters with wrong extent:", highcat(nrow(wrong_extent)), "/", highcat(nrow(sp_dirs_dt) * (length(hv.inc.t) + 1)))
      }
    } else {
      vebcat("Projection check completed successfully.", color = "proSuccess")
    }
  }

  catn("\nRunning statistics check.")

  # Check if stats has the same length as the number of species directories

  for (i in seq_along(stats_dir)) {
    stat_files <- list.files(stats_dir[[i]], full.names = TRUE)

    missing_sp_log <- paste0(logs_dir, "/", hv.method, "-missing-species.csv")

    dup_sp_log <- paste0(logs_dir, "/", hv.method, "-duplicate-species.csv")
    inc_dup_sp_log <- paste0(logs_dir, "/", hv.method, "-included-duplicate-species.csv")

    # For each stat file in the directory
    for (j in seq_along(stat_files)) {
      stat <- stat_files[[j]]

      vebcat("stat:", veb = verbose)
      vebprint(stat, veb = verbose)
      
      dt <- fread(stat)

      # Check for duplicates
      dt_dups <- dt[dt[, .I[duplicated(species)]]]
      fwrite(dt_dups, dup_sp_log, bom = TRUE)
      
      catn("Number of duplicated entries:", highcat(nrow(dt_dups)), "/", highcat(nrow(dt)))

      if (nrow(dt_dups) > 0) {
        vebcat("Removing", nrow(dt_dups), "duplicated species.", veb = verbose)
        dt <- dt[dt[, .I[!duplicated(species)]]]
      }

      dt <- dt[dt$excluded == FALSE, ]
      dt <- merge(dt, merged_dt[, c("species", "hvIteration")], by = "species", all.x = TRUE)
      dt_inc_dups <- dt[dt[, .I[duplicated(species)]]]

      if (nrow(dt) != nrow(merged_dt)) {
        vebcat(basename(stat), "failed the length check.", color = "nonFatalError", veb = verbose)
        
        catn("  Length of stats / length of projection dir:\n  ", highcat(nrow(dt)), "/", highcat(nrow(merged_dt)))
        missing_sp <- data.table(species = character(), dirName = character(), missingFrom = character(), hvIteration = numeric())

        # Create data tables for missing_sp_stat and missing_sp_merged if they are not empty
        missing_sp_stat <- dt[!dt$species %in% merged_dt$species, ]

        if (nrow(missing_sp_stat) > 0) {
          missing_sp_stat[, missingFrom := "projection"]
          missing_sp <- rbind(missing_sp, missing_sp_stat)
        }

        missing_sp_merged <- merged_dt[!merged_dt$species %in% dt$species, ]
        if (nrow(missing_sp_merged) > 0) {
          missing_sp_merged[, missingFrom := "stats"]
          missing_sp <- rbind(missing_sp, missing_sp_merged)
        }

        missing_sp <- merge(missing_sp, merged_dt[, c("species", "hvIteration")], by = "species", all.x = TRUE)
        missing_sp <- missing_sp[, -"hvIteration.y"]
        names(missing_sp) <- c("species", "dirName", "missingFrom", "hvIteration")
        
        vebcat("missing_sp:", veb = verbose)
        vebprint(missing_sp, veb = verbose)

        if (length(missing_sp) > 0) {
          vebcat(basename(stat), "Failed the length check.", color = "nonFatalError")
          fwrite(missing_sp, missing_sp_log, bom = TRUE)
        } else {
          vebcat(basename(stat), "passed the length check.", color = "proSuccess")
        }
      }
    }
  }
  
  vebcat("Hypervolume output check completed successfully", color = "funSuccess")

  return(list(
    missing_species = missing_sp,
    duplicate_species = dt_dups,
    missing_raster = missing_rast,
    wrong_crs = wrong_crs,
    wrong_extent = wrong_extent
  ))
}

clean_hv_output <- function(checked_data, out.dir, hv.dir, hv.projection, verbose = F) {
  vebcat("Initiating Hypervolume output cleaning", color = "funInit")

  # Check if the data is in the expected format
  if (length(checked_data) != 5) {
    stop("Data is not in the expected format. Expecting output from check_hv_output function.")

    for (i in seq_along(checked_data)) {
      data <- checked_data[[i]]

      if ("data.frame" %in% class(item) || "data.table" %in% class(item)) {
        vebcat("Found:", class(data), "expects a data.table or data.frame.", color = "fatalError")
        stop("Check your input data, expects output from the check_hv_output function")
      } else {
        stop("Data not a data.frame or data.table")
      }
    }
  } else {
    vebcat("Data is in the expected format.", veb = verbose)
  }

  # prepare logs
  logs_dir <- paste0(out.dir, "/logs/clean-data")
  create_dir_if(logs_dir)
  failed_reproj_log <- paste0(logs_dir, "/failed-projections.txt")
  create_file_if(failed_reproj_log)
  reproject_timer <- start_timer("reprojectRasters")

  # Clean the stats file
  stats_dir <- paste0(hv.dir, "/stats")

  stat_files <- list.files(stats_dir, full.names = TRUE)

  if (nrow(checked_data$duplicate_species) == 0) {
    catn("Stats file already clean.")
  } else {
    for (i in seq_along(stats_dir)) {
      stat_file_dir <- list.files(stats_dir[[i]], full.names = TRUE)
      stat_file <- stat_file_dir[[1]]

      stats_dt <- fread(stat_file)

      stats_dt_unique <- unique(stats_dt, by = "species", fromLast = TRUE)
      
      catn("Removing", highcat(nrow(stats_dt) - nrow(stats_dt_unique)), "species from the stats file and keeping the latest occurrence.")
    }
  }

  # clean extents
  if (nrow(checked_data$wrong_extent) == 0) {
    vebcat("extents already clean.", color = "proSuccess")
  } else {
    for (i in 1:nrow(checked_data$wrong_extent)) {
      rast_file <- checked_data$wrong_extent$fileName[i]

      cat("\rCropping", highcat(i), "/", highcat(nrow(checked_data$wrong_extent)))

      expected_res <- checked_data$wrong_extent$expectedRes[[i]]

      rast <- terra::rast(rast_file)

      extent <- checked_data$wrong_extent$expectedExt[[i]]

      extent <- strsplit(extent, ",")[[1]]

      extent <- as.numeric(extent)

      extent <- ext(extent)

      if (basename(rast_file) == "probability.tif") {
        new_ext <- terra::crop(rast, extent)
      } else {
        rast <- ceiling(rast)
        rast_new <- rast(res = expected_res, extent = extent)
        rast <- resample(rast, rast_new, method = "near")
        new_ext <- terra::crop(rast, ext(as.double(-180), as.double(180), as.double(55.7916666666667), as.double(83.6250000000)))
        ext(new_ext) <- extent
      }

      writeRaster(new_ext, rast_file, overwrite = TRUE)
    }
  }

  if (nrow(checked_data$wrong_crs) == 0) {
    catn("Projections already have the correct projections.")
  } else {
    for (i in 1:nrow(checked_data$wrong_crs)) {
      rast_file <- checked_data$wrong_crs$fileName[i]

      cat("\rReprojecting", highcat(i), "/", highcat(nrow(checked_data$wrong_crs)))

      rast <- terra::rast(rast_file)

      old_proj <- terra::crs(rast, proj = TRUE)

      if (!identical(old_proj, crs(hv.projection, proj = TRUE))) {
        if (basename(rast_file) == "probability.tif") {
          proj <- terra::project(rast, hv.projection, method = "bilinear")
        } else {
          proj <- terra::project(rast, hv.projection, method = "near")
        }

        new_proj <- terra::crs(proj, proj = TRUE)

        if (!identical(old_proj, new_proj)) {
          vebcat("Successfully reprojected", rast_file, color = "proSuccess", veb = verbose)
          writeRaster(proj, rast_file, overwrite = TRUE)
        } else {
          vebcat("Failed to reproject", rast_file, color = "nonFatalError", veb = verbose)
          try(frl <- file(failed_reproj_log, open = "at"))
          writeLines(failed_reproj_log, con = frl)
          close(frl)
        }
      } else {
        vebcat(rast_file, "Already in the correct projection.", color = "proSuccess", veb = verbose)
      }
    }
  }

  end_timer(reproject_timer)
  
  vebcat("Successfully completed Hypervolume output cleaning.", color = "funSuccess")

  return(stats_dt_unique)
}
