source_all("./src/hypervolume/data_analysis/components")

analyze_hv_stats <- function(region_hv, sp_hv, spec.name, verbose) {
  cat("Analyzing hypervolume for", cc$lightSteelBlue(sp_hv@Name), "species. \n")
  
  hv_set <- hypervolume_set(sp_hv, region_hv, check.memory = F)
  
  hv_stats <- hypervolume_overlap_statistics(hv_set)
  
  sp_surviv_region <- 1 - hv_stats[[4]]
  
  cat("Volume of", spec.name, "overlap in the CAVM", cc$lightSteelBlue(sp_surviv_region), "\n")

  return(hv_stats)
}

analyze_region_hv <- function(biovars, name, method, verbose) {
  cat(blue("Analyzing cavm hypervolume. \n"))

  directory <- paste0("./outputs/hypervolume/data_analysis/region/", tolower(name), "/")

  if (!dir.exists(directory)) dir.create(directory, recursive = T)

  if (!file.exists(paste0(directory, "hypervolume_", method, ".rds"))) {
    if (verbose) cat("File not found, initating hypervolume sequence. \n")
    matrix <- terra::values(biovars, na.rm = T)
    
    if (verbose) if (any(is.na(matrix))) cat(red("Some biovars region values are NA. \n")) else cat(green("No biovars region values are NA. \n"))

    if (verbose) cat("Matrix sample: \n")
    if (verbose) print(head(matrix, 3))

    hv <- hypervolume(matrix, name = name, method = method, verbose = verbose)

    # Save the hypervolume to a file
    saveRDS(hv, paste0(directory, "hypervolume_", method, ".rds"))
  } else {
    if (verbose) cat(name, "Hypervolume found, loading file. \n")
    # Load the hypervolume from the file
    hv <- readRDS(paste0(directory, "hypervolume_", method, ".rds"))
  }

  cat(cc$lightGreen(name, paste0("Hypervolume_", method), "Analysis Complete. \n"))
  return(hv)
}

plot_hypervolumes <- function(hv_list) {
  cat(blue("Plotting hypervolumes. \n"))
  for (i in seq_along(hv_list)) {
    hv <- hv_list[[i]]
    
    plot_dir <- "./outputs/data_analysis/hypervolume/species/plots"
    
    create_dir_if(plot_dir)
    
    png(paste0(plot_dir, "/hv_", names(hv_list)[i], ".png"), width = 800, height = 800, pointsize = 20)
    
    plot(hv, main = paste("Hypervolume for species", names(hv_list)[i]), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
    
    dev.off()
  }
}

plot_projects <- function(hv_list) {
  cat(blue("Plotting hypervolume_projections. \n"))

  for (i in seq_along(hv_list)) {
    hv <- hv_list[[i]]
    
    plot_dir <- "./outputs/data_analysis/hypervolume/species_projects/plots"
    create_dir_if(plot_dir)

    cat("Plotting", cc$lightSteelBlue(names(hv_list[[i]])), "\n")

    png(paste0(plot_dir, "/hv_", names(hv_list[[i]]), ".png"), width = 800, height = 800, pointsize = 20)

    plot(hv, main = paste("Hypervolume for", names(hv_list[[i]])), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)

    dev.off()
  }
}

get_ind_hv <- function(env_data) {
  cat(blue("Analyzing hypervolumes for"), cc$lightSteelBlue(length(env_data)), blue("species. \n"))
  
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

