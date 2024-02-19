visualize_histogram <- function(sp_cells, region, region.sub.color, region.sub.shades, region.name, plot.show = TRUE, verbose = T) {
  create_dir_if("./outputs/visualize/plots")
  
  # Figure 1A Whole CAVM
  cat(blue("Creating histogram for the entire Region\n"))
  
  print(class(sp_cells))
  print(head(sp_cells, 3))
  
  if (verbose) cat("Creating plot 1A. \n")
  fig1A <- ggplot(sp_cells, aes(x = cellSum) ) +
    # After_Stat accesses ggplot processed data, count represents the count of data points in each bin of the histogram
    geom_freqpoly(binwidth = 0.1, aes(y = after_stat(count)/sum(after_stat(count)))) +  
    labs(x = "Sum of Species Counts per Cell", y = "Proportion of Species Counts per Cell", title = paste0("Distribution of potential species richness in ", region.name)) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig1A)
  
  ggsave("./outputs/visualize/plots/figure-1A.png", fig1A, width = 10, height = 6)
  
  # Figure 1B different regions
  cat(blue("Creating histogram for each region in the", region.name, "\n"))
  
  # ADD shades of country for each floreg
  
  fig1B <- ggplot(sp_cells, aes(x = cellSum, color =  as.factor(sp_cells[[region.sub.color]]), linetype = as.factor(sp_cells[[region.sub.shades]])))  +
    geom_freqpoly(binwidth = 2, aes(y = after_stat(count)/sum(after_stat(count)))) +
    #scale_color_manual(values = unlist(color_list)) +
    labs(
      x = "Sum of Species Counts per Cell", y = "Proportion of Species Counts per Cell", 
      title = paste0("Distribution of potential species richness in different regions of ", region.name), 
      color = "Country", 
      linetype = "Floristic Province"
    ) +
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig1B)
  
  ggsave("./outputs/visualize/plots/figure-1B.png", fig1B, width = 10, height = 6)
  
  # Plot 1C - how many species [y] have a certain proportion of cells [x]
}

visualize_hotspots <- function(sp_cells, prob.value, region, region.sub, region.name, plot.show = TRUE, verbose = T) {

  #inc <- sp_cells[[1]]$sp_count
  #prob <- sp_cells[[2]]$sp_count

  #inc[inc == 0] <- NA
  #min_inc <- min(terra::values(inc$prop), na.rm = TRUE)
  #max_inc <- max(terra::values(inc$prop), na.rm = TRUE)
  min_inc <- min(sp_cells$cell_sum, na.rm = TRUE)
  max_inc <- max(sp_cells$cell_sum, na.rm = TRUE)
  region_ext <- terra::ext(region)
  
  fig2 <- ggplot() +
    geom_spatvector(data = region, aes_string(color = paste0("as.factor(", region.sub, ")"))) +
    scale_color_grey( guide = guide_legend(reverse = TRUE)) +
    #geom_spatraster(data = inc, aes(fill = prop)) +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, breaks = c(seq(min_inc, max_inc, by = 0.2)), limits = c(min_inc, max_inc), guide = guide_legend(reverse = TRUE)) +
    labs(x = "Longitude", y = "Latitude", title = paste0("Potential species distribution in the ", region.name), fill = "Proportion", color = "Country") +
    coord_sf(xlim = c(region_ext$xmin, region_ext$xmax), ylim = c(region_ext$ymin, region_ext$ymax)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig2)
  stop()
  
  ggsave("./outputs/visualize/plots/figure-2.png", plot = fig2)
  
  prob[prob == 0] <- NA
  prob_min <- min(terra::values(prob$proportion), na.rm = TRUE)
  prob_max <- max(terra::values(prob$proportion), na.rm = TRUE)
  fig3 <- ggplot() +
    geom_spatvector(data = cavm, aes_string(color = paste0("as.factor(", region.sub, ")")), fill = NA) +
    scale_color_grey( guide = guide_legend(reverse = TRUE)) +
    geom_spatraster(data = prob, aes_string(fill = prob.value)) +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, breaks = c(seq(prob_min, prob_max, by = 0.2)), limits = c(prob_min, prob_max), guide = guide_legend(reverse = TRUE)) +
    labs(x = "", y = "Latitude", title = paste0("Predicted species distribution in ", region.name), fill = "Relative Suitability Score", color = "Country") +
    coord_sf(xlim = c(region_ext$xmin, region_ext$xmax), ylim = c(region_ext$ymin, region_ext$ymax)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig3)
  
  ggsave("./outputs/visualize/plots/figure-3.png", plot = fig3)
}

visualize_matrix <- function(sp_cells, region.sub, plot.show = F) {
  cells <- sp_cells$region[!is.na(prop), ]
  
  min_prop <- min(cells$prop)
  max_prop <- max(cells$prop)
  
  break_scale <- seq(min_prop, 1, by = 0.10)
  cells$prop_bin <- cut(cells$prop, breaks = break_scale, labels = FALSE)
  
  fig4 <- ggplot(cells, aes(x = prop_bin, y = as.factor(cells[[region.sub]]), fill = prop)) +
    geom_tile(color = "white") +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, breaks = break_scale, guide = guide_legend(reverse = TRUE)) +
    labs(title = "Species climate suitability", x = "Species proportion", y = "Country") +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic")) +
    scale_x_discrete(breaks = break_scale)
    
  if (plot.show) print(fig4)

  ggsave("./outputs/visualize/plots/figure-4.png", plot = fig4)
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
  
  fig5 <- ggplot(data = sankey_dt, aes(axis1 = srcCountry, axis2 = floristicProvince, y = speciesSum)) +
    scale_x_discrete(limits = c("Origin", "Destination"), expand = c(.1, .1)) +
    xlab("Country") +
    ylab("Sums of each floristic region") +
    geom_alluvium(aes(fill = species)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 1.5) +
    theme_minimal() +
    ggtitle("Sankey plot from origin country to Arctic country")
  
  if (plot.show) print(fig5)
  
  ggsave("./outputs/visualize/plots/figure-5.png", device = "jpeg", unit = "px", width = 3840, height = 2160, plot = fig5)

  cat(cc$lightGreen("Data successfully visualized in a Sankey plot\n"))
}






