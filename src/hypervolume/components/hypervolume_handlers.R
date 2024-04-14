analyze_hv_stats <- function(region_hv, sp_hv, spec.name, verbose) {
  catn("Analyzing hypervolume for", highcat(sp_hv@Name), "species. \n")
  
  hv_set <- hypervolume_set(sp_hv, region_hv, check.memory = F)
  
  hv_stats <- hypervolume_overlap_statistics(hv_set)
  
  sp_surviv_region <- 1 - hv_stats[[4]]
  
  catn("Volume of", spec.name, "overlap in the CAVM", highcat(sp_surviv_region))
  
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
  }; catn()
  
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
  
  # Check column names
  missing_columns <- init.dt[!names(init.dt) %in% names(res)]
  extra_columns <- names(res)[!names(res) %in% names(init.dt)]
  
  if (extra_columns > 0) {
    catn("Length of output data table is higher than expected.")
    vebprint(extra_columns, text = "Extra columns:")
    check <- FALSE
  }
  
  if (missing_columns > 0) {
    catn("Length of output data table is lower than expected.")
    vebprint(missing_columns, text = "missing columns:")
    check <- FALSE
  }
  
  # Check duplicate entries
  stats_file <- paste0(hv.dir, "/stats.csv")
  
  existing_names <- unique(fread(stats_file, select = "cleanName"))
  
  out_name <- unique(res$cleanName)
  
  if (length(out_name) > 1) {
    catn("Multiple unique names found:", out_name)
    check <- FALSE
  }
  
  if (out_name %in% existing_names) {
    catn(out_name, "already in stats csv file.")
    check <- FALSE
  }
  
  # Check Projection outputs
  proj_dir <- paste0(hv.dir, "/projections/", out_name)
  all_files <- list.files(path = proj_dir, pattern = "\\.tif$", full.names = TRUE)
  prob_files <- list.files(path = proj_dir, pattern = "probability\\.tif$", full.names = TRUE)
  inc_files <- list.files(path = proj_dir, pattern = "inclusion\\.tif$", full.names = TRUE)
  
  # check dir length
  if (length(list.files(proj_dir)) > length(hv.incl.threshold) + 1) {
    catn("There are too many projection files in:", proj_dir)
    check <- FALSE
  } else if (length(list.files(proj_dir)) < length(hv.incl.threshold) + 1) {
    catn("There are too few projection files in:", proj_dir)
    check <- FALSE
  }
  
  # check crs
  for (i in 1:length(all_files)) {
    file <- all_files[[i]]
    
    r <- rast(file)
    
    r_crs <- terra::crs(r, proj = TRUE)
    
    e_crs <- terra::crs(longlat_crs, proj = TRUE)
    
    if (!identical(r_crs, e_crs)) {
      catn("The CRS is not identical:", proj_dir)
      vebprint(e_crs, text = "Expected CRS:")
      vebprint(r_crs, text = "output CRS:")
      check <- FALSE
    }
  }
  
  return(check)
} 