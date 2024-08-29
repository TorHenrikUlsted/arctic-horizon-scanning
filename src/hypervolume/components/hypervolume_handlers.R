analyze_hv_stats <- function(sp_hv, region_hv, spec.name, verbose) {
  catn("Analyzing hypervolume for", highcat(sp_hv@Name), "\n")

  # Check if the hypervolumes have the same number of dimensions
  if (length(sp_hv@Dimensionality) != length(region_hv@Dimensionality)) {
    stop("The species hypervolume and region hypervolume have different numbers of dimensions")
  }

  # Check if the hypervolumes have the same variable names
  if (!all(colnames(sp_hv@Data) %in% colnames(region_hv@Data))) {
    catn("Column names in species and region hypervolumes don't match.")
    vebprint(colnames(sp_hv@Data), text = "Species hypervolume names found:")
    vebprint(colnames(region_hv@Data), text = "Region hypervolume names found:")
    catn("Possible to rename with colnames(hypervolume@Data) <- new_names")
    stop("Column names in species and region hypervolumes don't match.")
  }

  if (!all(names(sp_hv@Parameters$kde.bandwidth) %in% names(region_hv@Parameters$kde.bandwidth))) {
    catn("Bandwidth names in species and region hypervolumes don't match.")
    vebprint(names(sp_hv@Parameters$kde.bandwidth), text = "Species Bandwidth names found:")
    vebprint(names(region_hv@Parameters$kde.bandwidth), text = "Region Bandwidth names found:")
    catn("Possible to rename with names(hypervolume@Parameters$kde.bandwidth) <- new_names")
    stop("Bandwidth names in species and region hypervolumes don't match.")
  }

  if (!all(colnames(sp_hv@RandomPoints) %in% colnames(region_hv@RandomPoints))) {
    catn("RandomPoint names in species and region hypervolumes don't match.")
    vebprint(colnames(sp_hv@RandomPoints), text = "Species RandomPoint names found:")
    vebprint(colnames(region_hv@RandomPoints), text = "Region RandomPoint names found:")
    catn("Possible to rename with names(hypervolume@Parameters$kde.bandwidth) <- new_names")
    stop("RandomPoint names in species and region hypervolumes don't match.")
  }

  catn("Using hypervolume_set")
  tryCatch(
    {
      hv_set <- hypervolume_set(sp_hv, region_hv, check.memory = FALSE)
    },
    error = function(e) {
      print(colnames(sp_hv@Data))
      print(colnames(region_hv@Data))
      print(colnames(hv_set@Data))
      catn("Error in hypervolume_set:", e$message)
      stop(e)
    }
  )

  catn("Using hypervolume_overlap_statistics")
  hv_stats <- hypervolume_overlap_statistics(hv_set)

  sp_surviv_region <- 1 - hv_stats[[4]]
  sp_realized_niche <- 1 - hv_stats[[3]]

  catn("Volume of", spec.name, "potential overlap in the region:", highcat(sp_surviv_region))
  catn("Volume of", spec.name, "Realized potential niche range:", highcat(sp_realized_niche))

  return(hv_stats)
}

analyze_region_hv <- function(biovars, file.out, method, verbose) {
  vebcat("Analyzing cavm hypervolume.", color = "funInit")

  create_dir_if(dirname(file.out))

  if (!file.exists(file.out)) {
    vebcat("File not found, analyzing hypervolume region.", veb = verbose)
    matrix <- terra::values(biovars, na.rm = TRUE)

    if (any(is.na(matrix))) {
      vebcat("Some biovars region values are NA.", color = "nonFatalError")
    } else {
      vebcat("No biovars region values are NA.", color = "proSuccess")
    }

    vebprint(head(matrix, 3), text = "Matrix sample:")

    hv <- hypervolume(
      matrix,
      name = "region",
      method = method,
      verbose = verbose
    )

    catn("Saving region hypervolume to:", colcat(file.out, color = "output"))

    saveRDS(hv, file.out)
  } else {
    catn("Hypervolume for region found, loading file.")
    # Load the hypervolume from the file
    hv <- readRDS(file.out)
  }

  vebcat(paste0("Hypervolume_", method), "Analysis Complete.", color = "funSuccess")

  return(hv)
}

analyze_ind_hv <- function(env_data, verbose = FALSE) {
  vebcat("Analyzing hypervolumes for", highcat(length(env_data)), "species.", color = "funInit")

  hv_list <- list()

  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(env_data), style = 3)

  for (species in names(env_data)) {
    catn("Creating hypervolume for species number", highcat(species))
    # Get the data for the species
    data <- env_data[[species]]

    vbprint(head(data, 3), veb = verbose)
    # Initialize an empty list to store the hypervolumes for the species
    species_hvs <- list()

    # For each dimension in the data
    for (dimension in colnames(data)) {
      catn("Using dimension:", highcat(dimension))
      print(head(data[, dimension]), 3)
      # Create a hypervolume for the dimension
      hv <- hypervolume_box(data[, dimension, drop = FALSE])
      catn("Add", dimension, "to species list.")
      # Add the hypervolume to the list
      species_hvs[[dimension]] <- hv

      vebcat("Adding dimension to list.", veb = verbose)
    }

    # Create a HypervolumeList for the species
    catn("Creating hypervolumeList")
    species_hv_list <- do.call(hypervolume_join, species_hvs)

    # Add the HypervolumeList to the list
    vebcat("Add species hypervolumeList to main list.", veb = verbose)
    hv_list[[species]] <- species_hv_list

    # Update the progress bar
    setTxtProgressBar(pb, which(names(env_data) == species))
  }
  catn()

  # Close the progress bar
  close(pb)

  vebcat("Analyzed hypervolumes for", highcat(length(env_data)), "species successfully", color = "funSuccess")

  return(hv_list)
}

plot_hypervolumes <- function(hv_list, out.dir) {
  vebcat("Plotting hypervolumes.", color = "funInit")

  for (i in seq_along(hv_list)) {
    hv <- hv_list[[i]]

    create_dir_if(out.dir)

    catn("Plotting", highcat(names(hv_list[[i]])))

    png(paste0(out.dir, "/hv_", names(hv_list)[i], ".png"), width = 800, height = 800, pointsize = 20)

    plot(hv, main = paste("Hypervolume for species", names(hv_list)[i]), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)

    dev.off()
  }

  vebcat("Plotting hypervolumes completed successfully", color = "funSuccess")
}

check_hv_results <- function(res, init.dt, hv.dir, hv.incl.threshold = 0.5, verbose = FALSE) {
  catn("Checking hypervolume output.")
  check <- TRUE
  vebprint(names(init.dt), verbose, "Init data table:")
  vebprint(names(res), verbose, "Res data table:")

  vebcat("Checking missing and extra columns.", veb = verbose)

  # Check column names
  missing_columns <- names(init.dt)[!names(init.dt) %in% names(res)]
  extra_columns <- names(res)[!names(res) %in% names(init.dt)]

  if (length(extra_columns) > 0 || length(missing_columns) > 0) {
    catn("Found", length(names(res)), "columns.")
    catn("Expected", length(names(init.dt)), "columns.")
  }

  if (length(extra_columns) > 0) {
    catn("Length of output data table is higher than expected.")
    vebprint(extra_columns, text = "Extra columns:")
    check <- FALSE
  }

  if (length(missing_columns) > 0) {
    catn("Length of output data table is lower than expected.")
    vebprint(missing_columns, text = "missing columns:")
    check <- FALSE
  }

  vebcat("Checking duplicate entires.", veb = verbose)

  # Check duplicate entries
  stats_file <- paste0(hv.dir, "/stats/stats.csv")

  existing_names <- unique(fread(stats_file, select = "cleanName")$cleanName)

  out_name <- unique(res$cleanName)

  vebprint(out_name, veb = verbose, text = "out_name")
  vebprint(existing_names, veb = verbose, text = "Existing names:")

  if (length(out_name) > 1) {
    catn("Multiple unique names found:", out_name)
    check <- FALSE
  }

  if (length(existing_names) != 0 & out_name %in% existing_names) {
    catn(out_name, "already in stats csv file.")
    check <- FALSE
  }

  vebcat("Finished checking duplicates.", veb = verbose)

  if (all(res$excluded == FALSE)) {
    vebcat("Checking projections.", veb = verbose)

    # Check Projection outputs
    proj_dir <- paste0(hv.dir, "/projections/", gsub(" ", config$species$file_separator, out_name))
    all_files <- list.files(path = proj_dir, pattern = "\\.tif$", full.names = TRUE)
    prob_files <- list.files(path = proj_dir, pattern = "probability\\.tif$", full.names = TRUE)
    inc_files <- list.files(path = proj_dir, pattern = "inclusion\\.tif$", full.names = TRUE)

    # check dir length
    if (length(all_files) > length(hv.incl.threshold) + 1 ||
      length(all_files) < length(hv.incl.threshold) + 1) {
      catn("Looked for files in:", proj_dir)
      catn("Found", length(all_files), "files.")
      catn("Expected", (length(hv.incl.threshold) + 1), "files.")
    }

    if (length(all_files) > length(hv.incl.threshold) + 1) {
      catn("There are too many projection files in:", proj_dir)
      check <- FALSE
    } else if (length(all_files) < length(hv.incl.threshold) + 1) {
      catn("There are too few projection files in:", proj_dir)
      check <- FALSE
    }

    # check crs
    for (i in 1:length(all_files)) {
      file <- all_files[[i]]

      print(file)

      r <- rast(file)

      print(r)

      r_crs <- terra::crs(r, proj = TRUE)

      e_crs <- terra::crs(get_crs_config(config$simulation$projection), proj = TRUE)

      if (!identical(r_crs, e_crs)) {
        catn("The CRS are not identical:", proj_dir)
        vebprint(e_crs, text = "Expected CRS:")
        vebprint(r_crs, text = "output CRS:")
        check <- FALSE
      }
    }
  } else if (!all(res$excluded == FALSE) & any(res$excluded == FALSE)) {
    catn("Some of the results, but not all, return FALSE. This is should not be possible..")
    check <- FALSE
  }

  if (check) catn("Finisihed hypervolume result check with NO faults.")
  if (!check) catn("Finisihed hypervolume result check WITH faults.")

  return(check)
}
