analyze_no_overlap_cases <- function(known.file, unknown.file, plot.dir = ".", vis.save.device = "jpeg", verbose = FALSE) {
  # Load both datasets
  known_stats <- fread(known.file)
  unknown_stats <- fread(unknown.file)
  
  # Get unique species
  known_unique <- unique(known_stats, by = "cleanName")
  unknown_unique <- unique(unknown_stats, by = "cleanName")
  
  # Function to extract and analyze no-overlap cases
  analyze_dataset <- function(stats_dt, dataset_name) {
    # Get species that failed due to no overlap
    no_overlap_species <- stats_dt[
      excluded == TRUE & 
        log(observations) > dimensions & 
        overlapRegion == 0
    ]
    
    # Basic statistics
    summary_stats <- no_overlap_species[, .(
      n_species = .N,
      mean_obs = mean(observations),
      median_obs = median(observations),
      min_obs = min(observations),
      max_obs = max(observations),
      mean_niche = mean(realizedNiche),
      mean_volume = mean(fracVolumeSpecies)
    )]
    
    # Taxonomic patterns
    taxon_patterns <- no_overlap_species[, .N, by = .(order, family)]
    setorder(taxon_patterns, -N)
    
    # Compare with successful species
    successful_species <- stats_dt[excluded == FALSE]
    
    wilcox_test <- wilcox.test(
      no_overlap_species$observations,
      successful_species$observations
    )
    
    # Add dataset identifier
    no_overlap_species[, dataset := dataset_name]
    
    return(list(
      species = no_overlap_species,
      summary = summary_stats,
      taxonomy = taxon_patterns,
      wilcox = wilcox_test
    ))
  }
  
  # Analyze both datasets
  known_analysis <- analyze_dataset(known_unique, "Arctic")
  unknown_analysis <- analyze_dataset(unknown_unique, "Non-Arctic")
  
  # Print analysis
  cli::cli_h2("Analysis of Species with No Hypervolume Overlap")
  
  # Arctic species results
  cli::cli_h3("Arctic Species")
  cli::cli_bullets(c(
    "*" = "Number of species: {.field {known_analysis$summary$n_species}}",
    "*" = "Observation statistics:",
    "  -" = "Mean: {.field {round(known_analysis$summary$mean_obs,1)}}",
    "  -" = "Median: {.field {round(known_analysis$summary$median_obs,1)}}",
    "  -" = "Range: {.field {round(known_analysis$summary$min_obs,1)}} - {.field {round(known_analysis$summary$max_obs,1)}}",
    "*" = "Wilcoxon test p-value (vs successful): {.field {format(known_analysis$wilcox$p.value, scientific = TRUE)}}"
  ))
  
  # Non-Arctic species results
  cli::cli_h3("Non-Arctic Species")
  cli::cli_bullets(c(
    "*" = "Number of species: {.field {unknown_analysis$summary$n_species}}",
    "*" = "Observation statistics:",
    "  -" = "Mean: {.field {round(unknown_analysis$summary$mean_obs,1)}}",
    "  -" = "Median: {.field {round(unknown_analysis$summary$median_obs,1)}}",
    "  -" = "Range: {.field {round(unknown_analysis$summary$min_obs,1)}} - {.field {round(unknown_analysis$summary$max_obs,1)}}",
    "*" = "Wilcoxon test p-value (vs successful): {.field {format(unknown_analysis$wilcox$p.value, scientific = TRUE)}}"
  ))
  
  # Combine data for visualization
  plot_data <- rbindlist(list(
    known_analysis$species,
    unknown_analysis$species
  ))
  
  # Create visualizations
  p1 <- ggplot(plot_data) +
    geom_boxplot(aes(x = dataset, y = log(observations))) +
    labs(
      title = "Distribution of Observations for Species with No Overlap",
      x = "Dataset",
      y = "Log(Number of Observations)"
    ) +
    theme_minimal() +
    ggplot.theme()
  
  p2 <- ggplot(plot_data) +
    geom_point(aes(x = observations, y = realizedNiche, color = dataset)) +
    scale_x_log10() +
    labs(
      title = "Relationship between Observations and Realized Niche",
      x = "Number of Observations (log scale)",
      y = "Realized Niche Volume"
    ) +
    theme_minimal() +
    ggplot.theme()
  
  # Save plots
  save_ggplot(
    p1, 
    "no_overlap_observations", 
    3000, 3000, 
    plot.dir, 
    vis.save.device,
    suppress = TRUE
  )
  
  save_ggplot(
    p2, 
    "no_overlap_niche_relationship", 
    3000, 3000, 
    plot.dir, 
    vis.save.device,
    suppress = TRUE
  )
  
  return(list(
    arctic = known_analysis,
    non_arctic = unknown_analysis,
    plots = list(p1, p2)
  ))
}

analyze_patterns <- function(known.file, unknown.file, plot.dir = ".", vis.save.device = "jpeg", verbose = FALSE) {
  # Load and prepare data
  known_stats <- fread(known.file)
  unknown_stats <- fread(unknown.file)
  
  # Get unique species and identify no-overlap cases
  analyze_dataset <- function(stats_dt, dataset_name) {
    unique_dt <- unique(stats_dt, by = "cleanName")
    
    cli::cli_h3(paste("Analyzing", dataset_name))
    cli::cli_progress_bar(
      name = paste(dataset_name, "analysis"),
      type = "tasks",
      total = 8
    )
    
    # Group creation
    cli::cli_progress_update(
      status = "Creating analysis groups",
      set = 1
    )
    unique_dt[, group := fifelse(
      excluded == TRUE & log(observations) > dimensions & overlapRegion == 0,
      "No Overlap",
      "Successful"
    )]
    
    cli::cli_progress_update(
      status = "Separating cases",
      set = 2
    )
    # Separate into no-overlap and successful cases
    no_overlap <- unique_dt[
      excluded == TRUE & 
        log(observations) > dimensions & 
        overlapRegion == 0
    ]
    successful <- unique_dt[excluded == FALSE]
    
    # Data combination
    cli::cli_progress_update(
      status = "Preparing analysis dataset",
      set = 3
    )
    analysis_dt <- rbind(no_overlap, successful)
    analysis_dt[, group := fifelse(cleanName %in% no_overlap$cleanName, "No Overlap", "Successful")]
    
    # Geographic matrix
    cli::cli_progress_update(
      status = "Computing geographic distances",
      set = 4
    )
    geo_matrix <- as.matrix(analysis_dt[, .(level3Long, level3Lat)])
    geo_dist <- vegan::vegdist(geo_matrix, method = "euclidean")
    
    # Geographic testing
    cli::cli_progress_update(
      status = "Testing geographic patterns",
      set = 5
    )
    geo_test <- vegan::adonis2(
      geo_dist ~ group,
      data = analysis_dt,
      permutations = 999
    )
    
    # Geographic summary
    cli::cli_progress_update(
      status = "Summarizing geographic patterns",
      set = 6
    )
    geo_summary <- rbind(no_overlap, successful)[, .(
      mean_long = mean(level3Long),
      mean_lat = mean(level3Lat),
      sd_long = sd(level3Long),
      sd_lat = sd(level3Lat)
    ), by = group]
    
    wm <- get_world_map("longlat")
    
    
    # Plot geographic patterns
    p1 <- ggplot() +
      geom_spatvector(data = wm) +
      geom_point(data = rbind(no_overlap, successful), 
                 aes(x = level3Long, y = level3Lat, color = group), size = 0.75) +
      geom_point(data = geo_summary,
                 aes(x = mean_long, y = mean_lat, fill = group),
                 size = 4, shape = 23) +
      labs(
        title = paste(dataset_name, "- Geographic Distribution"),
        x = "Longitude",
        y = "Latitude"
      ) +
      theme_minimal() +
      ggplot.theme()
    
    # Save plots
    save_ggplot(
      p1, 
      paste0(tolower(dataset_name), "_geographic_patterns"), 
      3000, 3000, 
      save.dir = plot.dir, 
      save.device = vis.save.device,
      suppress = TRUE
    )
    
    # Environmental analysis
    cli::cli_progress_update(
      status = "Analyzing environmental patterns",
      set = 7
    )
    env_vars <- grep("^bio", names(unique_dt), value = TRUE)
    if (length(env_vars) == 0) {
      cli_alert_warning("No bioclimatic variables found in the dataset")
      env_test <- NULL
      p2 <- NULL
    } else {
      env_matrix <- as.matrix(analysis_dt[, ..env_vars])
      env_dist <- vegan::vegdist(env_matrix, method = "euclidean")
      
      env_test <- vegan::adonis2(
        env_dist ~ group,
        data = analysis_dt,
        permutations = 999
      )
      
      # Create environmental plot only if we have variables
      env_data <- melt(unique_dt,
                       id.vars = c("cleanName", "excluded"),
                       measure.vars = env_vars)
      
      p2 <- ggplot(env_data, aes(x = variable, y = value, fill = as.factor(excluded))) +
        geom_boxplot() +
        coord_flip() +
        labs(
          title = paste(dataset_name, "- Environmental Variables"),
          x = "Bioclimatic Variable",
          y = "Value",
          fill = "Excluded"
        ) +
        theme_minimal() +
        ggplot.theme()
     
     save_ggplot(
       p2, 
       paste0(tolower(dataset_name), "_environmental_patterns"),
       3000, 3000, 
       save.dir = plot.dir, 
       save.device = vis.save.device,
       suppress = TRUE
     )
   }
    
    cli::cli_progress_done()
    
    return(list(
      geographic_test = geo_test,
      environmental_test = env_test,
      geographic_summary = geo_summary,
      plots = list(geographic = p1, environmental = p2)
    ))
  }
  
  # Analyze both datasets
  arctic_patterns <- analyze_dataset(known_stats, "Arctic")
  non_arctic_patterns <- analyze_dataset(unknown_stats, "Non-Arctic")
  
  # Print results
  cli::cli_h2("Geographic and Environmental Pattern Analysis")
  
  # In the print section of analyze_patterns(), add null checks:
  
  cli_h3("Arctic Species")
  cli_bullets(c(
    "*" = "Geographic patterns (PERMANOVA):",
    "  -" = sprintf("F-statistic: %.3f", arctic_patterns$geographic_test$F[1]),
    "  -" = sprintf("P-value: %.2e", arctic_patterns$geographic_test$`Pr(>F)`[1]), 
    "  -" = sprintf("R-squared: %.3f", arctic_patterns$geographic_test$R2[1])
  ))
  
  if (!is.null(arctic_patterns$environmental_test)) {
    cli_bullets(c(
      "*" = "Environmental patterns (PERMANOVA):",
      "  -" = sprintf("F-statistic: %.3f", arctic_patterns$environmental_test$F[1]),
      "  -" = sprintf("P-value: %.2e", arctic_patterns$environmental_test$`Pr(>F)`[1]),
      "  -" = sprintf("R-squared: %.3f", arctic_patterns$environmental_test$R2[1])
    ))
  }
  
  cli_h3("Non-Arctic Species") 
  cli_bullets(c(
    "*" = "Geographic patterns (PERMANOVA):",
    "  -" = sprintf("F-statistic: %.3f", non_arctic_patterns$geographic_test$F[1]),
    "  -" = sprintf("P-value: %.2e", non_arctic_patterns$geographic_test$`Pr(>F)`[1]),
    "  -" = sprintf("R-squared: %.3f", non_arctic_patterns$geographic_test$R2[1])
  ))
  
  if (!is.null(non_arctic_patterns$environmental_test)) {
    cli_bullets(c(
      "*" = "Environmental patterns (PERMANOVA):",
      "  -" = sprintf("F-statistic: %.3f", non_arctic_patterns$environmental_test$F[1]),
      "  -" = sprintf("P-value: %.2e", non_arctic_patterns$environmental_test$`Pr(>F)`[1]),
      "  -" = sprintf("R-squared: %.3f", non_arctic_patterns$environmental_test$R2[1])
    ))
  }
  
  return(list(
    arctic = arctic_patterns,
    non_arctic = non_arctic_patterns
  ))
}
