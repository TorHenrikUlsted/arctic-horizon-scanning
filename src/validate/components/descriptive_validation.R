find_dataset_name <- function(dataset = "known", compare, example = FALSE, verbose = FALSE) {
  
  spec_known <- config$dataset$known
  spec_unknown <- config$dataset$unknown
  
  if (dataset == "unknown") {
    compare[, sourceDataset := spec_unknown]
    return(compare)
  }
  
  res_dir <- file.path("./outputs/validation", spec_unknown, "results")
  create_dir_if(res_dir)
  
  file <- paste0(
    "./outputs/filter/", 
    if (dataset == "known") spec_known else spec_unknown, 
    "/filtered-sp-keys.csv"
  )
  
  data <- fread(file)
  
  res <- unique(compare, by = "scientificName")
  
  catn("input res_keys:", highcat(nrow(res)))
  
  res_keys <- get_gbif_keys(
    spec = res$scientificName,
    out.dir = file.path(res_dir, "hv-data")
  )
  
  res_keys <- res_keys[, .(cleanName, scientificName, speciesKey, usageKey)]
  
  res <- res_keys[res, on = "scientificName"]
  
  catn("output res_keys:", highcat(nrow(res_keys)))
  
  setnames(data, c("speciesKey.GBIF"), c("speciesKey"))
  
  data <- data[, .(speciesKey, usageKey, scientificName, species, sourceDataset)]
  
  # Add cleanName
  data[, cleanName := {
    res <- lapply(cli_progress_along(scientificName, "downloading"), function(i) {
      clean_spec_name(scientificName[i], config$species$standard_symbols, config$species$standard_infraEpithets, verbose)
    })
    catn()
    vapply(res, function(x) x$cleanName, character(1))
  }]
  
  catn("Initial length:", highcat(nrow(res)))
  
  matched_key <- data[res, on = "usageKey", nomatch = 0L]
  catn("Matched after usageKey:", highcat(nrow(matched_key)))
  
  unmatched_key <- res[!(usageKey %in% matched_key$usageKey)]
  catn("Unmatched after usageKey:", highcat(nrow(unmatched_key)))
  
  # Filter data for unmatched rows
  data_unmatched <- data[!(usageKey %in% matched_key$usageKey)]
  catn("Data left after usageKey:", highcat(nrow(data_unmatched)))
  
  
  
  matched_cleanName <- data_unmatched[unmatched_key, on = "cleanName", nomatch = 0L]
  catn("Matched after cleanName:", highcat(nrow(matched_cleanName)))
  
  unmatched_cleanName <- unmatched_key[!(cleanName %in% matched_cleanName$cleanName)]
  catn("Unmatched after cleanName:", highcat(nrow(unmatched_cleanName)))
  
  # Filter data for unmatched rows
  data_unmatched <- data_unmatched[!(cleanName %in% matched_cleanName$cleanName)]
  catn("Data left after cleanName:", highcat(nrow(data_unmatched)))
  
  unmatched_cleanName[, `:=` (
    verbatimCleanName = cleanName,
    cleanName = i.cleanName,
    i.cleanName = NULL
  )]
  
  matched_icleanName <- data_unmatched[unmatched_cleanName, on = "cleanName", nomatch = 0L]
  catn("Matched after icleanName:", highcat(nrow(matched_icleanName)))
  
  unmatched_icleanName <- unmatched_cleanName[!(cleanName %in% matched_icleanName$cleanName)]
  catn("Unmatched after icleanName:", highcat(nrow(unmatched_icleanName)))
  
  # Filter data for unmatched rows
  data_unmatched <- data_unmatched[!(cleanName %in% matched_icleanName$cleanName)]
  catn("Data left after icleanName:", highcat(nrow(data_unmatched)))
  
  
  matched_spKey <- unique(data_unmatched[unmatched_icleanName, on = "speciesKey", nomatch = 0L], by = "scientificName")
  matched_spKey <- unique(matched_spKey, by = "verbatimCleanName")
  catn("Matched after speciesKey:", highcat(nrow(matched_spKey)))
  
  unmatched_spKey <- unmatched_icleanName[!(speciesKey %in% matched_spKey$speciesKey)]
  catn("Unmatched after speciesKey:", highcat(nrow(unmatched_spKey)))
  
  all_matches <- rbind(
    matched_key,
    matched_cleanName,
    matched_icleanName,
    matched_spKey,
    fill = TRUE
  )
  
  all_matches <- unique(all_matches, by = "scientificName")

  
  # Find columns starting with "i."
  i_cols <- grep("^i\\.", names(all_matches), value = TRUE)
  
  # Remove these columns
  all_matches[, (i_cols) := NULL]
  
  catn("All results combined:", highcat(nrow(all_matches)))
  dsets <- unique(all_matches$sourceDataset)
  catn("Found", highcat(length(dsets)), "unique datasets")
  if (length(dsets) < 20) {
    catn()
    vebprint(dsets, text = "Datasets:")
  }
  
  return(all_matches)
}

add_too_few_observations <- function(dir, stats, verbose = FALSE) {
  spec <- fread(file.path(dir, "removed-species.csv"))
  
  if (grepl("_validation", dir)) {
    combined <- rbind(stats, spec, fill = TRUE)
  } else {
    setnames(spec, "species", "cleanName")
    
    spec[, excluded := TRUE]
    
    combined <- rbind(stats, spec, fill = TRUE)
  }
  
  return(combined)
}

analyze_hypervolume_exclusions <- function(known.dir, unknown.dir, stats.dir, example = FALSE, plot.dir = ".", vis.save.device = "jpeg", verbose = FALSE) {
  cli::cli_alert_info("Loading datasets")
  known_stats <- fread(file.path(known.dir, "stats.csv"))
  unknown_stats <- fread(file.path(unknown.dir, "stats.csv"))
  
  known_stats <- add_too_few_observations(
    dir = known.dir, 
    stats = known_stats, 
    verbose = verbose
  )
  
  known_unique <- find_dataset_name(
    dataset = "known",
    compare = known_stats,
    example = example
  )
  
  unknown_stats <- add_too_few_observations(
    dir = unknown.dir,
    stats = unknown_stats,
    verbose = verbose
  )
  
  unknown_unique <- find_dataset_name(
    dataset = "unknown",
    compare = unknown_stats,
    example = example
  )
  
  unknown_unique <- unique(unknown_unique, by = "cleanName")
  
  catn(
    sum(known_unique$excluded == TRUE), 
    "/", sum(known_unique$excluded == FALSE),
    "/", nrow(known_unique)
  )
  catn(
    sum(unknown_unique$excluded == TRUE), 
    "/", sum(unknown_unique$excluded == FALSE), 
    "/", nrow(unknown_unique)
  )

  cli::cli_h2("Dataset Overview")
  cli::cli_alert_info("Checking duplicates due to multiple botanical regions")
  cli::cli_bullets(c(
    "*" = "Arctic dataset: 
    {.field {nrow(known_unique)}} records reduced to 
    {.field {nrow(known_unique[excluded == FALSE])}} unique species",
    "*" = "Non-Arctic dataset: 
    {.field {nrow(unknown_unique)}} records reduced to 
    {.field {nrow(unknown_unique[excluded == FALSE])}} unique species"
  ))

  # Function to analyze exclusions for a dataset
  analyze_dataset <- function(stats_dt, dataset_name) {
    cols_to_select <- c(
      "cleanName",
      "scientificName",
      "observations",
      "dimensions",
      "excluded",
      "jaccard",
      "overlapRegion",
      "realizedNiche",
      "sourceDataset"
    )
    
    # Initialize results tracking
    results <- stats_dt[, ..cols_to_select]

    results[, logObservations := log(observations)]


    # Add exclusion reason
    results[, exclusionReason := fifelse(
      !excluded, "Included",
      fifelse(
        logObservations <= dimensions,
        "Too Few Observations",
        fifelse(
          overlapRegion == 0,
          "No Hypervolume Overlap",
          "Other"
        )
      )
    )]


    # Get excluded species details
    excluded_species <- results[excluded == TRUE]
    
    calc_stats <- function(dataset.name, dataset) {
      # Calculate summary stats
      total <- nrow(dataset)
      excluded <- sum(dataset$excluded)
      included <- sum(!dataset$excluded)
      too_few <- sum(dataset$exclusionReason == "Too Few Observations")
      no_overlap <- sum(dataset$exclusionReason == "No Hypervolume Overlap")
      
      cli::cli_h3(dataset.name)
      cli::cli_bullets(c(
        "*" = "Total species analyzed: {.field {total}}",
        
        "*" = "Total included: {.field {included}
      ({round(100 * included/total, 1)}%)}",
      
      "*" = "Total excluded: {.field {excluded}
      ({round(100 * excluded/total, 1)}%)}",
      
      "*" = "Excluded due to too few observations: 
      {.field {too_few} ({round(100 * too_few/total, 1)}%)}",
      
      "*" = "Excluded due to no hypervolume overlap: 
      {.field {no_overlap} ({round(100 * no_overlap/total, 1)}%)}"
      ))
    }
    
    calc_stats(
      dataset.name = dataset_name,
      dataset = results
    )
    
    # G
    sub_datasets <- unique(results$sourceDataset)

    if (length(sub_datasets) > 1) {
      for (i in 1:length(sub_datasets)) {
        name <- sub_datasets[i]
        dt <- results[sourceDataset == name]
        calc_stats(
          dataset.name = toTitleCase(name),
          dataset = dt
        )
      }
    }

    return(list(
      results = results,
      species.excluded = excluded_species
    ))
  }
  
 

  cli::cli_h2("Exclusion Analysis")
  known_analysis <- analyze_dataset(known_unique, "Arctic Species")
  unknown_analysis <- analyze_dataset(unknown_unique, "Non-Arctic Species")
  
  # Generate visualization if any species were excluded
  if (nrow(known_analysis$species.excluded) > 0 || nrow(unknown_analysis$species.excluded)) {
    # Combine data for plotting
    known_analysis$species.excluded[, dataset := "Arctic"]
    unknown_analysis$species.excluded[, dataset := "Non-Arctic"]
    plot_data <- rbindlist(list(
      known_analysis$species.excluded,
      unknown_analysis$species.excluded
    ))

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
  
  known_analysis$results[, dataset := "Arctic"]
  unknown_analysis$results[, dataset := "Non-Arctic"]
  
  plot_data <- rbindlist(list(
    known_analysis$results,
    unknown_analysis$results
  ))
  

  plot_data[, proportion := {
    reason <- exclusionReason 
    fifelse(
      reason == "Included",
      sum(reason == "Included", na.rm = TRUE) / .N,
      fifelse(
        reason == "Too Few Observations",
        sum(reason == "Too Few Observations", na.rm = TRUE) / .N,
        fifelse(
          reason == "No Hypervolume Overlap",
          sum(reason == "No Hypervolume Overlap", na.rm = TRUE) / .N,
          NA
        )
      )
    )
  }, by = .(dataset, sourceDataset)]
  
  plot_data[, N := .N, by = .(dataset, sourceDataset)]
  
  plot_summary <- plot_data[, .(proportion = first(proportion)), 
                            by = .(sourceDataset, exclusionReason, N)]
  
  # fwrite(
  #   plot_data[exclusionReason == "No Hypervolume Overlap" & sourceDataset == "aba", .(cleanName, scientificName)],
  #   "./aba-species.csv",
  #   bom = TRUE)
  
  # Create exclusion reason plot
  p3 <- ggplot(plot_summary, aes(x = sourceDataset, y = proportion, fill = exclusionReason)) +
    geom_bar(stat = "identity", position = "fill") +
    geom_text(
      aes(x = sourceDataset, y = 1.05, label = paste0("N = ", N)),
      inherit.aes = FALSE,
      size = 3
    ) +
    ggplot.filler(
      gradient = config$ggplot$gradient$vis.gradient,
      begin = 0.94,
      end = 0,
      guide = list(ncol = 1, nrow = 3),
      scale.type = "fill-d",
      na.value = config$ggplot$gradient$na.value
    ) +
    labs(
      title = "",
      x = "Dataset",
      y = "Proportion",
      fill = "Exclusion Reason"
    ) +
    theme_minimal() +
    ggplot.theme()
  
  save_ggplot(
    p3,
    "spec-exclusion-reason",
    3000, 3000,
    save.dir = plot.dir,
    save.device = vis.save.device,
    suppress = TRUE
  )

  # Return results invisibly
  invisible(list(
    known = known_analysis,
    unknown = unknown_analysis
  ))
}

analyze_hypervolume_validation <- function(known.dir, unknown.dir, example = FALSE, plot.dir = ".", vis.save.device = "jpeg", verbose = FALSE) {
  cli::cli_h1("Hypervolume Validation Analysis")

  # First run the exclusions analysis
  exclusions <- analyze_hypervolume_exclusions(
    known.dir = known.dir,
    unknown.dir = unknown.dir,
    example = example,
    plot.dir = plot.dir,
    vis.save.device = vis.save.device,
    verbose = FALSE
  )
  
  known_results <- exclusions$known$results
  unknown_results <- exclusions$unknown$results

  # Calculate model validation metrics
  known_res <- nrow(known_results)
  correctly_identified_known <- sum(!known_results$excluded)
  known_few <- sum(known_results$exclusionReason == "Too Few Observations")
  missed_known <- sum(known_results$exclusionReason == "No Hypervolume Overlap")
  
  unknown_res <- nrow(unknown_results)
  potential_new_known <- sum(!unknown_results$excluded)
  unknown_few <- sum(unknown_results$exclusionReason == "Too Few Observations")
  missed_unknown <- sum(unknown_results$exclusionReason == "No Hypervolume Overlap")
  
  cli::cli_h2("Model Validation Results")

  # Known Arctic Species Performance
  cli::cli_h3("Performance on Known Arctic Species")
  region_detection_rate <- correctly_identified_known / (known_res - known_few) * 100
  
  cli::cli_bullets(c(
    "*" = "Total known Arctic species: {.field {known_res}}",
    
    "*" = "Correctly identified: 
    {.field {correctly_identified_known}}",
    
    "*" = "Failed to identify: {.field {missed_known}} species",
    
    "*" = "Success rate of analyzed species: {.field ({round(region_detection_rate, 1)}%)}",
    
    "*" = "This suggests the model is 
    {.field {ifelse(region_detection_rate > 80, 
    'is reliable', 'may need adjustment')}} 
    at identifying Arctic-suitable climates"
  ))

  # Potential New Arctic Species
  cli::cli_h3("Potential New Arctic Species Analysis")
  cli::cli_bullets(c(
    "*" = "Identified {.field {potential_new_known}} 
    non-Arctic species with climate niche overlap",
    
    "*" = "Success rate of analyzed species: {.field 
    {round(potential_new_known/(unknown_res - unknown_few)*100, 1)}%}",
    
    "*" = "These species merit further investigation for invasion potential"
  ))

  if (!is.null(plot.dir)) {
    create_dir_if(plot.dir)

    # Create scatter plot of niche volume vs region overlap
    p2 <- ggplot() +
      geom_point(
        data = known_results,
        aes(x = realizedNiche, y = overlapRegion, color = "Arctic")
      ) +
      geom_point(
        data = unknown_results,
        aes(x = realizedNiche, y = overlapRegion, color = "Non-Arctic")
      ) +
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
  
  
  ## Plot stacked barplot
  p5 <- ggplot() +
    geom_bar(data = known_results, aes(x = overlapRegion, fill = "Arctic"), alpha = 0.5) +

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
      percentage = potential_new_known / nrow(unknown_results) * 100
    )
  ))
}

