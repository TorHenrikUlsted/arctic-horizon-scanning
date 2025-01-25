analyze_hypervolume_exclusions <- function(known.file, unknown.file, plot.dir = ".", vis.save.device = "jpeg", verbose = FALSE) {
  
  cli::cli_alert_info("Loading datasets")
  known_stats <- fread(known.file)
  unknown_stats <- fread(unknown.file)
  
  # Get unique species
  known_unique <- unique(known_stats, by = "cleanName")
  unknown_unique <- unique(unknown_stats, by = "cleanName")
  
  cli::cli_h2("Dataset Overview")
  cli::cli_alert_info("Checking duplicates due to multiple botanical regions")
  cli::cli_bullets(c(
    "*" = "Arctic dataset: {.field {nrow(known_stats)}} records reduced to {.field {nrow(known_unique)}} unique species",
    "*" = "Non-Arctic dataset: {.field {nrow(unknown_stats)}} records reduced to {.field {nrow(unknown_unique)}} unique species"
  ))
  
  # Function to analyze exclusions for a dataset
  analyze_dataset <- function(stats_dt, dataset_name) {
    cols_to_select <- c(
      "cleanName",
      "observations",
      "dimensions",
      "excluded",
      "jaccard",
      "overlapRegion",
      "realizedNiche"
    )
    
    # Initialize results tracking
    results <- stats_dt[, ..cols_to_select]
    
    results[, logObservations := log(observations)]

    
    # Add exclusion reason
    results[, exclusionReason := fifelse(
      !excluded, "Not Excluded",
      fifelse(logObservations <= dimensions, 
              "Too Few Observations",
              fifelse(overlapRegion == 0, 
                      "No Hypervolume Overlap",
                      "Other"))
    )]
    
    
    # Get excluded species details
    excluded_species <- results[excluded == TRUE]
    
    # Calculate summary stats
    total <- nrow(results)
    excluded <- sum(results$excluded)
    too_few <- sum(results$exclusionReason == "Too Few Observations")
    no_overlap <- sum(results$exclusionReason == "No Hypervolume Overlap")
    
    cli::cli_h3(dataset_name)
    cli::cli_bullets(c(
      "*" = "Total species analyzed: {.field {total}}",
      "*" = "Total excluded: {.field {excluded}} ({round(100 * excluded/total, 1)}%)",
      "*" = "Excluded due to too few observations: {.field {too_few} ({round(100 * too_few/total, 1)}%)}",
      "*" = "Excluded due to no hypervolume overlap: {.field {no_overlap} ({round(100 * no_overlap/total, 1)}%)}"
    ))
    
    return(list(
      results = results,
      species.excluded = excluded_species
    ))
  }
  
  cli::cli_h2("Exclusion Analysis")
  known_analysis <- analyze_dataset(known_unique, "Arctic Species")
  unknown_analysis <- analyze_dataset(unknown_unique, "Non-Arctic Species")
  
  # Generate visualization if any species were excluded
  if(nrow(known_analysis$species.excluded) > 0 || nrow(unknown_analysis$species.excluded)) {
    # Combine data for plotting
    known_analysis$species.excluded[, dataset := "Arctic"]
    unknown_analysis$species.excluded[, dataset := "Non-Arctic"]
    plot_data <- rbindlist(list(
      known_analysis$species.excluded,
      unknown_analysis$species.excluded
    ))
    
    # Create exclusion reason plot
    p <- ggplot(plot_data, aes(x = dataset, fill = exclusionReason)) +
      geom_bar(position = "dodge") +
      scale_fill_manual(values = c("Too Few Observations" = "#4682B4", "No Hypervolume Overlap" = "#b4464e")) +
      labs(
        title = "Reasons for Species Exclusion",
        x = "Dataset",
        y = "Number of Species",
        fill = "Exclusion Reason"
      ) +
      theme_minimal() +
      ggplot.theme()
    
    save_ggplot(
      p,
      "spec-exclusion-reason",
      3000, 3000,
      save.dir = plot.dir,
      save.device = vis.save.device,
      suppress = TRUE
    )
    
    # Create observations vs dimensions scatter plot
    p2 <- ggplot(plot_data, aes(x = logObservations, y = dimensions, color = exclusionReason)) +
      geom_point() +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      labs(
        title = "Log(Observations) vs Dimensions for Excluded Species",
        x = "Log(Observations)",
        y = "Number of Dimensions",
        color = "Exclusion Reason"
      ) +
      theme_minimal() +
      ggplot.theme()
    
    save_ggplot(
      p2,
      "nobs_dimensions",
      3000, 3000,
      save.dir = plot.dir,
      save.device = vis.save.device,
      suppress = TRUE
    )
  }
  
  # Return results invisibly
  invisible(list(
    known = known_analysis,
    unknown = unknown_analysis
  ))
}


analyze_hypervolume_validation <- function(known.file, unknown.file, plot.dir = ".", vis.save.device = "jpeg", verbose = FALSE) {
  cli::cli_h1("Hypervolume Validation Analysis")
  
  # First run the exclusions analysis
  exclusions <- analyze_hypervolume_exclusions(
    known.file = known.file,
    unknown.file = unknown.file,
    plot.dir = plot.dir,
    vis.save.device = vis.save.device,
    verbose = FALSE
  )
  
  known_results <- exclusions$known$results
  unknown_results <- exclusions$unknown$results
  
  # Calculate model validation metrics
  known_res <- nrow(known_results)
  correctly_identified_known <- sum(!known_results$excluded)
  missed_known <- sum(known_results$excluded)
  potential_new_known <- sum(!unknown_results$excluded)
  
  cli::cli_h2("Model Validation Results")
  
  # Known Arctic Species Performance
  cli::cli_h3("Performance on Known Arctic Species")
  region_detection_rate <- correctly_identified_known / known_res * 100
  cli::cli_bullets(c(
    "*" = "Total known Arctic species: {.field {known_res}}",
    "*" = "Correctly identified: {.field {correctly_identified_known}} ({round(region_detection_rate, 1)}%)",
    "*" = "Failed to identify: {.field {missed_known}} species",
    "*" = "This suggests the model is {.field {ifelse(region_detection_rate > 80, 'reliable', 'may need adjustment')}} at identifying Arctic-suitable climates"
  ))
  
  # Analysis of exclusion reasons for missed Arctic species
  if(missed_known > 0) {
    missed <- known_results[excluded == TRUE]
    cli::cli_h3("Analysis of Missed Arctic Species")
    missed_reasons <- missed[, .N, by = exclusionReason]
    cli::cli_bullets(c(
      "*" = missed_reasons[, paste0(exclusionReason, ": {.field ", N, "}")] 
    ))
  }
  
  # Potential New Arctic Species
  cli::cli_h3("Potential New Arctic Species Analysis")
  cli::cli_bullets(c(
    "*" = "Identified {.field {potential_new_known}} non-Arctic species with climate niche overlap",
    "*" = "This represents {.field {round(potential_new_known/nrow(unknown_results)*100, 1)}%} of analyzed non-Arctic species",
    "*" = "These species merit further investigation for invasion potential"
  ))
  
  if (!is.null(plot.dir)) {
    create_dir_if(plot.dir)
    
    # Create scatter plot of niche volume vs region overlap
    p2 <- ggplot() +
      geom_point(data = known_results, 
                 aes(x = realizedNiche, y = overlapRegion, color = "Arctic")) +
      geom_point(data = unknown_results, 
                 aes(x = realizedNiche, y = overlapRegion, color = "Non-Arctic")) +
      labs(
        title = "Realized Niche Volume vs Region Overlap",
        x = "Realized Niche Volume",
        y = "Region Overlap"
      ) +
      scale_color_manual(values = c("Arctic" = "#4682B4", "Non-Arctic" = "#b4464e")) +
      theme_minimal() +
      ggplot.theme()
    
      save_ggplot(
        p2,
        "volume_overlap_relationship",
        3000, 3000,
        save.dir = plot.dir,
        save.device = vis.save.device,
        suppress = TRUE
      )
    
    # Create visualization comparing niche volumes
    p3 <- ggplot() +
      geom_density(data = known_results, aes(x = realizedNiche, fill = "Arctic"), alpha = 0.5) +
      geom_density(data = unknown_results, aes(x = realizedNiche, fill = "Non-Arctic"), alpha = 0.5) +
      labs(
        title = "Distribution of Realized Niche Volumes",
        x = "Realized Niche Volume",
        y = "Density"
      ) +
      scale_fill_manual(values = c("Arctic" = "#4682B4", "Non-Arctic" = "#b4464e")) +
      theme_minimal() +
      ggplot.theme()
    
    save_ggplot(
      p3,
      "realized_niche_volume_comparison",
      3000, 3000,
      save.dir = plot.dir,
      save.device = vis.save.device,
      suppress = TRUE
    )
    
    # Create visualization comparing niche volumes
    p4 <- ggplot() +
      geom_density(data = known_results, aes(x = overlapRegion, fill = "Arctic"), alpha = 0.5) +
      geom_density(data = unknown_results, aes(x = overlapRegion, fill = "Non-Arctic"), alpha = 0.5) +
      labs(
        title = "Distribution of Overlap Volumes",
        x = "Overlap Volume",
        y = "Density"
      ) +
      scale_fill_manual(values = c("Arctic" = "#4682B4", "Non-Arctic" = "#b4464e")) +
      theme_minimal() +
      ggplot.theme()
    
    save_ggplot(
      p4,
      "overlap_volume_comparison",
      3000, 3000,
      save.dir = plot.dir,
      save.device = vis.save.device,
      suppress = TRUE
    )
  }
  
  # Return results invisibly
  invisible(list(
    known_performance = list(
      total = known_res,
      identified = correctly_identified_known,
      missed = missed_known,
      detection.rate = region_detection_rate
    ),
    potential_new = list(
      count = potential_new_known,
      percentage = potential_new_known/nrow(unknown_results)*100
    )
  ))
}