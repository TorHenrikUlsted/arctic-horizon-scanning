visualize_freqpoly <- function(sp_cells, region, region.name, plot.x, plot.y, plot.color, plot.shade, plot.show = FALSE, verbose = FALSE) {
  create_dir_if("./outputs/visualize/plots")
  
  # Figure 1A Whole CAVM
  cat(blue("Creating histogram for the entire Region\n"))
  
  print(class(sp_cells))
  print(head(sp_cells, 3))
  
  if (verbose) cat("Creating plot 1A. \n")
  
  sp_cells[, richness := ifelse(richness == 0, NA, richness)]
  
  fig1A <- ggplot(sp_cells, aes_string(x = plot.x) ) +
    # After_Stat accesses ggplot processed data, count represents the count of data points in each bin of the histogram
    geom_freqpoly(binwidth = 0.1, aes(y = ..count.. / sum(..count..)) ) +
    labs(
      x = "Species Richness (log10)", 
      y = "Relative Cell Frequency", 
      title = paste0("Potential Species Richness in ", region.name)
      ) +
    scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig1A)
  
  ggsave("./outputs/visualize/plots/figure-1A.jpeg", device = "jpeg", unit = "px", width = 2160, height = 2160, fig1A)
  
  # Figure 1B different regions
  cat(blue("Creating histogram for each region in the", region.name, "\n"))
  
  # Order by plot.color
  sp_cells <- sp_cells[order(sp_cells[[plot.color]]), ]
  sp_cells[[plot.color]] <- factor(sp_cells[[plot.color]], levels = unique(sp_cells[[plot.color]]))
  
  # ADD shades of country for each floreg
  plcol <- eval(as.factor(sp_cells[[plot.color]]))
  plshd <- eval(as.factor(sp_cells[[plot.shade]]))
  
  fig1B <- ggplot(sp_cells, aes_string(x = plot.x, color = plcol, fill = plshd))  +
    # ..count.. / sum(..count..) calculates the proportion of cells at each species richness value for each region within each bin of the histogram
    geom_freqpoly(binwidth = 0.1, aes(y =  ..count.. / sum(..count..)) ) + 
    scale_color_viridis_d(guide = "legend", option = "B") + 
    labs(
      x = "Species Richness (log10)", 
      y = "Relative Cell Frequency", 
      title = paste0("Potential Species Richness in Different Regions of ", region.name), 
      color = "Country", 
      linetype = "Floristic Province"
    ) +
    scale_y_continuous(limits = c(0, NA), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig1B)
  
  ggsave("./outputs/visualize/plots/figure-1B.jpeg",device = "jpeg", unit = "px", width = 2160, height = 2160, fig1B)
  
  # Plot 1C - how many species [y] have a certain proportion of cells [x]
}

visualize_hotspots <- function(rast, region, region.name, extent, projection, projection.method, plot.show = FALSE, verbose = FALSE) {
  cat(blue("Visualizing Potential Species hotspots\n"))
  
  cat("Checking crs.\n")
  rast <- check_crs(rast, projection, projection.method)
  region <- check_crs(region, projection, projection.method)

  cat("Getting min and max values.\n")
  min_lim <- where.min(rast)
  min_lim <- min(min_lim[, 3])
  max_lim <- where.max(rast)
  max_lim <- max(max_lim[, 3])
  
  cat("Plotting hotspots.\n")
  
  fig2A <- ggplot() +
    geom_spatvector(data = world_map) +
    geom_spatraster(data = rast) +
    scale_fill_whitebox_c(
      palette = "muted", 
      guide = guide_legend(reverse = TRUE), 
      limits = (c(min_lim, max_lim)), 
      #trans = "log",
      labels = function(x) format(x, big.mark = ",", scientific = FALSE)
    ) +
    labs(title = paste0("Potential species hotspots in the CAVM"), fill = "Species Richness") +
    coord_sf(xlim = c(extent$xmin, extent$xmax), ylim = c(extent$ymin, extent$ymax)) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 12, face = "bold.italic"),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 8)
    )
  
  if (plot.show) print(fig2A)
  
  cat("Saving plot.\n")
  ggsave("./outputs/visualize/plots/figure-2A.jpeg", device = "jpeg", unit = "px", width = 2160, height = 2160, plot = fig2A)
  
  cat(cc$lightGreen("Successfully visualized Potential Species hotspots\n"))
}

visualize_highest_spread <- function(rast, region, region.name, extent, projection, projection.method, plot.show = FALSE, verbose = FALSE) {
  cat(blue("Visualizing Potential Species hotspots\n"))
  
  cat("Checking crs.\n")
  rast <- check_crs(rast, projection, projection.method, verbose = verbose)
  region <- check_crs(region, projection, projection.method, verbose = verbose)  
  
  cat("Plotting hotspots.\n")
  
  fig2B <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = rast) +
    facet_wrap(~lyr, ncol = 3, labeller = label_wrap_gen(width = 20)) +
    scale_fill_whitebox_c(palette = "muted", breaks = c(0, 1), labels = c("No Overlap", "Overlap"), guide = guide_legend(reverse = TRUE)) +
    labs(title = paste0("Species with Highest Potential Spread in the CAVM"), fill = "Overlap Value") +
    coord_sf(xlim = c(extent$xmin, extent$xmax), ylim = c(extent$ymin, extent$ymax)) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, size = 12, face = "bold.italic"),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 8)
    )
  
  if (plot.show) print(fig2B)
  
  cat("Saving plot.\n")
  ggsave("./outputs/visualize/plots/figure-2B.jpeg", device = "jpeg", unit = "px",  plot = fig2B)
  
  cat(cc$lightGreen("Successfully visualized Potential Species hotspots\n"))
}

visualize_suitability <- function(rast, region, region.name, plot.show = F, verbose = F) {
  
  if (!identical(crs(rast, proj = TRUE), crs(region, proj = TRUE))) {
    cat("Reprojecting to laea.\n")
    cat(crs(rast, proj = TRUE), "\n")
    cat(crs(region, proj = TRUE), "\n")
    cat(identical(crs(rast, proj = TRUE), crs(region, proj = TRUE)), "\n")
    rast <- project(rast, laea_crs, method = "bilinear")
  }
  
  cat("Acquiring min and max values.\n")
  prob_min <- 0.001
  prob_max <- where.max(rast[[1]])[[3]]
  region_ext <- ext(region)
  
  cat("Generating plot.\n")
  fig3 <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = rast) +
    scale_fill_whitebox_c(palette = "muted", limits = c(prob_min, prob_max), breaks = c(seq(100, prob_max, by = 100)), guide = guide_legend(reverse = TRUE)) +
    facet_wrap(~lyr, nrow = 3) +
    labs(x = "Longitude", y = "Latitude", title = paste0("Predicted species distribution in ", region.name), fill = "Suitability Score", color = "Country") +
    coord_sf(xlim = c(region_ext$xmin, region_ext$xmax), ylim = c(region_ext$ymin, region_ext$ymax)) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 10, face = "bold.italic"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 10)
      )
  
  if (plot.show) print(fig3)
  ggsave("./outputs/visualize/plots/figure-3.jpeg", device = "jpeg", unit = "px", width = 2160, height = 2160, plot = fig3)
}

calculate_richness <- function(dt, sp_cols, taxon, verbose = F) {
  
  cat("Calculatig potential richness.\n")
  # Reform it to be richness and not abundance
  dt <- dt[, (sp_cols) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = sp_cols]
  
  # Calculate the sum for each FLOREG and keep specific columns
  dt <- dt[, c(lapply(.SD, sum), .(country = country, floristicProvince = floristicProvince)), by = FLOREG, .SDcols = sp_cols]
  
  dt <- unique(dt, by = "FLOREG")
  
  dt <- dt[, (sp_cols) := lapply(.SD, function(x) ifelse(x > 0, 1, x)), .SDcols = sp_cols]
  
  dt[, speciesRichness := rowSums(.SD > 0), .SDcols = sp_cols]
  
  dt <- melt(dt, id.vars = c("FLOREG", "country", "floristicProvince", "speciesRichness"), measure.vars = sp_cols, variable.name = "species", value.name = "count")
  print(dt)
  # Get higher taxon ranks
  cat("Getting higher taxon names.\n")
  checklist <- name_backbone_checklist(name_data = levels(dt$species))
  
  checklist <- data.table(checklist[, c("kingdom", "phylum", "order", "family", "genus", "canonicalName", "verbatim_name")])
  
  names(checklist)[7] <- "species"
  
  na_sp <- checklist[is.na(checklist$genus),]
  
  cat("Species missing classification.\n")
  print(na_sp)
  
  cat("Merging data tables.\n")
  merged_dt <- merge(dt, checklist, by = "species")
  
  dt <- data.table(merged_dt)
  
  dt <- dt[complete.cases(dt),]
  
  taxon_col <- dt[[taxon]]
  
  #dt <- dt[, .(country, floristicProvince, potentialRichness = length(unique(species))), by = .(FLOREG, taxon_col)]

  dt <- unique(dt, by = c(taxon, "FLOREG"))
  
  dt <- dt[, .(country, floristicProvince, kingdom, phylum, order, family, speciesRichness, totalRichness = sum(speciesRichness)), by = FLOREG]
  
  cat("Calculating relative richness\n")
  dt[, relativeRichness := speciesRichness / totalRichness]
  
  return(dt)
}

visualize_richness <- function(dt, axis.x, axis.y, fill, group, plot.show = F, verbose = F) {
  # Reorder the bars
  dt[[axis.x]] <- as.factor(dt[[axis.x]])
  
  if (verbose) {
    cat("Before sorting:\n")
    print(levels(dt[[axis.x]]))
  }
  
  # Order the entire data framy by the group parameter
  dt <- dt[order(dt[[group]]), ]
  # Order the levels of the x axis by the sorted dt 
  dt[[axis.x]] <- factor(dt[[axis.x]], levels = unique(dt[[axis.x]]))
  
  if (verbose) {
    cat("After sorting:\n")
    print(levels(dt[[axis.x]]))
  }
  
  fig4 <- ggplot(dt, aes(x = dt[[axis.x]], y = dt[[axis.y]], fill = dt[[fill]] )) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_whitebox_d(palette = "muted") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Floristic Region", y = "Relative Richness", fill = paste0(toupper(substr(fill, 1, 1)), substr(fill, 2, nchar(fill))))
  
  if (plot.show) print(fig4)
  
  ggsave("./outputs/visualize/plots/figure-4.jpeg", plot = fig4, device = "jpeg", unit="px")
  
}

visualize_sankey <- function(dt.src, dt.target, plot.show = F, verbose = F) {
  cat(blue("Visualizing data in a sankey plot\n"))
  # dt check
  if (!("data.table" %in% class(dt.src)) || !("data.table" %in% class(dt.target))) {
    stop("input data is not a data.table.")
  }
  
  cat("Setting up target data.\n")
  
  if (verbose) cat("Reshaping target data.\n")
  
  # Reshape to long format for sankey plotting
  target_dt <- melt(dt.target, id.vars = c("ID", "FLOREG", "country", "floristicProvince"), measure.vars = sp_cols, variable.name = "species", value.name = "count")
  # Sum species values
  target_dt <- target_dt[, speciesSum := sum(count, na.rm = TRUE), by = FLOREG]
  
  target_dt <- unique(target_dt, by = c("species", "FLOREG"))
  
  # Remove speciesSum of 0
  target_dt <- target_dt[speciesSum > 0]
  
  cat("Setting up source data.\n")
  # Get source country
  source_dt <- visualize_data$included_sp %>% 
    select(species, srcCountryIso = countryIso, srcCountry = country)
  
  # Calculate source country
  source_dt <- source_dt[, srcSpeciesSum := uniqueN(species), by = .(srcCountryIso, srcCountry)]
  
  source_dt <- unique(source_dt, by = c("srcCountryIso", "srcCountry"))
  
  source_dt <- source_dt[srcSpeciesSum > 0]
  
  # Merge the dts
  sankey_dt <- merge(target_dt, source_dt, by = "species", allow.cartesian = TRUE)
  
  ################
  #     Plot     #
  ################
  
  cat("Creating sankey plot.\n")
  
  fig5 <- ggplot(data = sankey_dt, aes(axis1 = srcCountry, axis2 = country, y = speciesSum)) +
    scale_x_discrete(limits = c("Origin", "Destination"), expand = c(.1, .1)) +
    xlab("Country") +
    ylab("Sums of each floristic region") +
    geom_alluvium(aes(fill = floristicProvince)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.5) +
    theme_minimal() +
    ggtitle("Sankey plot from origin country to Arctic country")
  
  if (plot.show) print(fig5)
  
  ggsave("./outputs/visualize/plots/figure-5.jpeg", device = "jpeg", unit = "px", width = 3840, height = 2160, plot = fig5)

  cat(cc$lightGreen("Data successfully visualized in a Sankey plot\n"))
}






