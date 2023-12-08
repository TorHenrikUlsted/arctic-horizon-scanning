source_all("./src/hypervolume/data_analysis/components")

analyze_correlation <- function(wc_region, threshold) {
  cat(blue("Analyzing correlation data. \n"))

  cor_mat <- get_correlation(wc_region, threshold)

  # Use the function
  # imp_biovars <- important_biovars(cor_mat, threshold)

  cat(cc$lightGreen("Correlation analysis completed. \n"))

  # return(imp_biovars)
}

analyze_hypervolume <- function(region_hv, sp_hv, verbose) {
  if (verbose == T) cat("Analyzing hypervolume for", cc$lightSteelBlue(sp_hv@Name), "species. \n")
  
  hv_set <- hypervolume_set(sp_hv, region_hv, check.memory = F)
  
  hv_stats <- hypervolume_overlap_statistics(hv_set)
  
  sp_surviv_region <- 1 - hv_stats[[4]]
  
  cat("Volume of species overlap in the CAVM", cc$lightSteelBlue(sp_surviv_region), "\n")

  return(sp_surviv_region)
}

analyze_region_hv <- function(biovars, name, method, samples.per.point, verbose) {
  cat(blue("Analyzing cavm hypervolume. \n"))

  directory <- paste0("./outputs/data_analysis/hypervolume/region/", tolower(name), "/")

  if (!dir.exists(directory)) dir.create(directory, recursive = T)

  if (!file.exists(paste0(directory, "hypervolume_", method, ".rds"))) {
    cat("File not found, initating hypervolume sequence. \n")
    matrix <- terra::values(biovars, na.rm = T)
    
    if (any(is.na(matrix))) cat(red("Some biovars region values are NA. \n")) else cat(green("No biovars region values are NA. \n"))

    cat("Matrix sample: \n")
    print(head(matrix, 3))

    hv <- hypervolume(matrix, name = name, method = method, samples.per.point = samples.per.point, verbose = verbose)

    # Save the hypervolume to a file
    saveRDS(hv, paste0(directory, "hypervolume_", method, ".rds"))
  } else {
    cat(name, "Hypervolume found, loading file. \n")
    # Load the hypervolume from the file
    hv <- readRDS(paste0(directory, "hypervolume_", method, ".rds"))
  }

  cat(cc$lightGreen(name, paste0("Hypervolume_", method), "Analysis Complete. \n"))
  return(hv)
}


get_ind_hv <- function(env_data) {
  cat("Analyzing hypervolumes for", cc$lightSteelBlue(length(env_data)), "species. \n")

  hv_list <- list()

  # Create a progress bar
  pb <- txtProgressBar(min = 0, max = length(env_data), style = 3)

  for (species in names(env_data)) {
    cat("Creating hypervolume for species number", cc$lightSteelBlue(species), "\n")
    # Get the data for the species
    data <- env_data[[species]]

    print(head(data, 3))
    # Initialize an empty list to store the hypervolumes for the species
    species_hvs <- list()

    # For each dimension in the data
    for (dimension in colnames(data)) {
      cat("Using dimension:", cc$lightSteelBlue(dimension), "\n")
      print(head(data[, dimension]), 3)
      # Create a hypervolume for the dimension
      hv <- hypervolume_box(data[, dimension, drop = FALSE])
      cat("Add", dimension, "to species list. \n")
      # Add the hypervolume to the list
      species_hvs[[dimension]] <- hv

      cat("Adding dimension to list \n")
    }

    # Create a HypervolumeList for the species
    cat("Creating hypervolumeList \n")
    species_hv_list <- do.call(hypervolume_join, species_hvs)

    # Add the HypervolumeList to the list
    cat("Add species hypervolumeList to main list \n")
    hv_list[[species]] <- species_hv_list

    # Update the progress bar
    setTxtProgressBar(pb, which(names(env_data) == species))
    cat("\n")
  }

  # Close the progress bar
  close(pb)

  return(hv_list)
}

analyze_inclusion = function(region_hv, sp_data, verbose) {
  if (verbose == T) cat("Running inclusion test for", cc$lightSteelBlue(names(sp_data)[1])) 
  if (verbose == T) cat("Samples per point", cc$lightSteelBlue(ceiling((10^(3 + sqrt(ncol(sp_data[[1]]))))/nrow(sp_data[[1]]))), "\n")
  
  inc_test <- hypervolume_inclusion_test(region_hv, 
                                         sp_data, 
                                         reduction.factor = 1, 
                                         fast.or.accurate = "accurate", 
                                         fast.method.distance.factor = 1,
                                         accurate.method.threshold = quantile(region_hv@ValueAtRandomPoints, 0.5), 
                                         verbose = verbose
                                         )
  if (verbose == T) cat(cc$lightGreen("inclusion analysis completed successfully \n"))
  
  return(inc_test)
}

plot_hypervolumes <- function(hv_list) {
  cat("Plotting hypervolumes. \n")
  for (i in seq_along(hv_list)) {
    hv <- hv_list[[i]]
    
    if (!dir.exists("./outputs/data_analysis/hypervolume/species/plots/")) dir.create("./outputs/data_analysis/hypervolume/species/plots/", recursive = T)
    
    png(paste0("./outputs/data_analysis/hypervolume/species/plots/hv_", names(hv_list)[i], ".png"), width = 800, height = 800, pointsize = 20)
    
    plot(hv, main = paste("Hypervolume for species", names(hv_list)[i]), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
    
    dev.off()
  }
}

plot_projects <- function(hv_list) {
  cat("Plotting hypervolume_projections. \n")

  for (i in seq_along(hv_list)) {
    hv <- hv_list[[i]]

    cat("Plotting", cc$lightSteelBlue(names(hv_list[[i]])), "\n")

    if (!dir.exists("./outputs/data_analysis/hypervolume/species_projects/plots/")) dir.create("./outputs/data_analysis/hypervolume/species_projects/plots/", recursive = T)

    png(paste0("./outputs/data_analysis/hypervolume/species_projects/plots/hv_", names(hv_list[[i]]), ".png"), width = 800, height = 800, pointsize = 20)

    plot(hv, main = paste("Hypervolume for Saxifraga oppositifolia", names(hv_list[[i]])), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)

    dev.off()
  }
}
