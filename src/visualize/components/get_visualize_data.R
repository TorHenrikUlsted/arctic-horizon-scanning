get_visualize_data <- function(hv.dir, hv.method, verbose = F, warn, err) {
  cat(blue("Initializing visulasation data. \n"))
  # Define directories
  vis_stats_dir <- "./outputs/visualize/stats"

  sp_dirs <- list.dirs(paste0(hv.dir, "/projections/", hv.method))
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]

  if (verbose) {
    cat("Sample directory names: \n")
    print(head(sp_dirs, 3))
  }

  if (!file.exists("./outputs/visualize/stats/included-species.csv") || !file.exists("./outputs/visualize/stats/included-species.csv")) {
    stats_src <- paste0(hv.dir, "/stats/box-stats.csv")
    sp_stats <- fread(stats_src)

    # Get the region for the species
    glonaf_wfo_one <- fread("./resources/synonym-checked/glonaf-species-wfo-one.csv")

    gwo <- glonaf_wfo_one %>%
      select(rawName.ORIG) %>%
      # Get the refined scientificNames as used in the hypervolume method
      mutate(species = apply(glonaf_wfo_one[, c("genus", "specificEpithet", "infraspecificEpithet"), drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " ")))

    gwo$species <- trimws(gwo$species)
    # Swap space with line to match the sp-stats
    gwo$species <- gsub(" ", "-", gwo$species)

    colnames(gwo) <- c("origName", "species")

    # Remove duplicates
    gwo <- gwo[!duplicated(gwo$species)]

    # Match the names with the sp_stats
    gwo_stats <- sp_stats %>%
      left_join(gwo, by = "species")

    # Get the region id
    glonaf_orig_df <- fread("./resources/data-raw/glonaf.csv")

    godf <- glonaf_orig_df %>%
      select(standardized_name, region_id)

    colnames(godf) <- c("origName", "regionId")

    glonaf_orig_region <- fread("./resources/region/glonaf-region-desc.csv")

    gor <- glonaf_orig_region %>% select(region_id, country_ISO, country)

    colnames(gor) <- c("regionId", "countryIso", "country")

    # Match the gwo_stats and godf --- Some will have multiple region ids
    matched_glonaf <- merge(gwo_stats, godf, by = "origName")

    if (verbose) cat("Number of rows when appending region ids:", cc$lightSteelBlue(nrow(matched_glonaf)), "\n")

    # Remove duplicates based on combined region id and standardized name
    matched_dups <- nrow(matched_glonaf[duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ])
    if (verbose) cat("Number of duplicates after appending region id:", cc$lightSteelBlue(matched_dups), "\n")

    if (matched_dups > 0) {
      matched_glonaf <- matched_glonaf[!duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ]

      if (verbose) cat("Number of duplicates after removal:", cc$lightSteelBlue(nrow(matched_glonaf[duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ])))

      if (verbose) cat("Number of rows after finishing appending region ids:", cc$lightSteelBlue(nrow(matched_glonaf)), "\n")
    }

    # Append region info
    matched_stats <- merge(matched_glonaf, gor, by = "regionId")

    matched_stats <- matched_stats %>%
      select(species, origName, observations, dimensions, samplesPerPoint, randomPoints, excluded, jaccard, sorensen, fracVolumeSpecies, fracVolumeRegion, overlapRegion, includedOverlap, regionId, countryIso, country)

    # Seperate into included and excluded species

    included_sp <- matched_stats %>% filter(excluded == FALSE)

    if (verbose) cat("Number of included species:", cc$lightSteelBlue(nrow(included_sp)), "\n")

    create_dir_if("./outputs/visualize/stats")
    fwrite(included_sp, paste0(vis_stats_dir, "/included-species.csv"), bom = TRUE)

    excluded_sp <- matched_stats %>% filter(excluded == TRUE)

    if (verbose) cat("Number of included species:", cc$lightSteelBlue(nrow(excluded_sp)), "\n")

    fwrite(excluded_sp, paste0(vis_stats_dir, "/excluded-species.csv"), bom = TRUE)
  } else {
    included_sp <- fread(paste0(vis_stats_dir, "/included-species.csv"))

    excluded_sp <- fread(paste0(vis_stats_dir, "/excluded-species.csv"))
  }

  # Get tifs
  inc_stack <- terra::rast()
  prob_stack <- terra::rast()
  
  withCallingHandlers(
    {
      for (i in seq_along(sp_dirs)) {
        species <- sp_dirs[[i]]

        cat("\rStacking raster files ", i, " / ", length(sp_dirs))

        flush.console()

        sp_inc <- paste0(species, "/inclusion-0.5.tif")
        sp_prob <- paste0(species, "/probability.tif")

        try(
          {
            sp_rast <- terra::rast(sp_inc)
            names(sp_rast) <- basename(species)
            inc_stack <- c(inc_stack, sp_rast)
            
            sp_rast <- terra::rast(sp_prob)
            names(sp_rast) <- basename(species)
            prob_stack <- c(prob_stack, sp_rast)
          },
          silent = TRUE
        )
      }
    },
    warning = function(w) warn(w, warn_txt = "Warning when stacking inclusion analyses."),
    error = function(e) err(e, err_txt = "Error when stacking inclusion analyses.")
  )

  cat(cc$lightGreen("\nvisualization data successfully initialised.\n"))

  return(list(
    included_sp = included_sp,
    excluded_sp = excluded_sp,
    inc_stack = inc_stack,
    prob_stack = prob_stack
  ))
}
