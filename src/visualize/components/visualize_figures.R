visualize_freqpoly <- function(sp_cells, region, region.name, plot.x, plot.y, plot.color, plot.shade, plot.show = FALSE, verbose = FALSE) {
  create_dir_if("./outputs/visualize/plots")
  
  # Figure 1A Whole CAVM
  vebcat("Creating histogram for the entire Region", color = "funInit")
  
  print(class(sp_cells))
  print(head(sp_cells, 3))
  
  vebcat("Creating plot 1A.", veb = verbose)
  
  #sp_cells[, richness := ifelse(richness == 0, NA, richness)]
  # Remove NA values
  sp_cells <- sp_cells[!is.na(sp_cells[[plot.x]]), ]
  
  # Remove Inf values
  sp_cells <- sp_cells[sp_cells[[plot.x]] != Inf, ]

  
  fig1A <- ggplot(sp_cells, aes_string(x = plot.x) ) +
    # After_Stat accesses ggplot processed data, count represents the count of data points in each bin of the histogram
    geom_freqpoly(binwidth = 0.1, aes(y = after_stat(count) / sum(after_stat(count))) ) +
    geom_text() +
    labs(
      x = "Species Richness (log10)", 
      y = "Proportion of cells in Arctic CAVM", 
      title = paste0("Potential Door-knocker Species Richness in ", region.name)
      ) +
    scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig1A)
  
  ggsave("./outputs/visualize/plots/figure-1A.jpeg", device = "jpeg", unit = "px", width = 2160, height = 2160, fig1A)
  
  vebcat("Frequency for the entire region successfully visualized", color = "funSuccess")
  
  # Figure 1B different regions
  vebcat("Creating histogram for each region in the", region.name, color = "funInit")
  
  # Order by plot.color
  sp_cells <- sp_cells[order(sp_cells[[plot.color]]), ]
  sp_cells[[plot.color]] <- factor(sp_cells[[plot.color]], levels = unique(sp_cells[[plot.color]]))
  
  # ADD shades of country for each floreg
  plcol <- eval(as.factor(sp_cells[[plot.color]]))
  plshd <- eval(as.factor(sp_cells[[plot.shade]]))
  
  fig1B <- ggplot(sp_cells, aes_string(x = plot.x, color = plcol, fill = plshd))  +
    # ..count.. / sum(..count..) calculates the proportion of cells at each species richness value for each region within each bin of the histogram
    geom_freqpoly(binwidth = 0.1, aes(y =  after_stat(count) / sum(after_stat(count))) ) + 
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
  
  vebcat("Frequency for each floristic region successfully visualized", color = "funSuccess")
}

visualize_hotspots <- function(rast, region, region.name, extent, projection, projection.method, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Potential Species hotspots", color = "funInit")
  
  catn("Checking crs.")
  rast <- check_crs(rast, projection, projection.method)
  region <- check_crs(region, projection, projection.method)

  catn("Getting min and max values.")
  min_lim <- where.min(rast)
  min_lim <- min(min_lim[, 3])
  max_lim <- where.max(rast)
  max_lim <- max(max_lim[, 3])
  catn("Plotting hotspots.")
  
  fig2A <- ggplot() +
    geom_spatvector(data = world_map) +
    geom_spatraster(data = rast) +
    scale_fill_viridis_b(
      option = "B", 
      #palette = "E",
      guide = guide_legend(reverse = TRUE),
      #direction = -1,
      #trans = "log1p",
      limits = (c(min_lim, max_lim)),
      breaks = c(0,1,5,10,50,100,500,1000, max_lim),
      labels = function(x) format((x), big.mark = ",", scientific = FALSE, digits = 2),
      na.value = "transparent"
    ) +
    labs(title = paste0("Potential door-knocker species hotspots in the CAVM"), fill = "Species Richness") +
    coord_sf(xlim = c(extent$xmin, extent$xmax), ylim = c(extent$ymin, extent$ymax)) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 12, face = "bold.italic"),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 8)
    )
  
  if (plot.show) print(fig2A)
  
  catn("Saving plot.")
  ggsave("./outputs/visualize/plots/figure-2A.jpeg", device = "jpeg", unit = "px", width = 2160, height = 2160, plot = fig2A)
  
  vebcat("Successfully visualized Potential Species hotspots", color = "funSuccess")
}

visualize_highest_spread <- function(rast, region, region.name, extent, projection, projection.method, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing highest spread", color = "funInit")
  
  catn("Checking crs.")
  rast <- check_crs(rast, projection, projection.method, verbose = verbose)
  region <- check_crs(region, projection, projection.method, verbose = verbose)  
  
  catn("Plotting hotspots.")
  
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
  
  fig2B_out <- "./outputs/visualize/plots/figure-2B.jpeg"
  
  catn("Saving plot to:", colcat(fig2B_out, color = "output"))
  ggsave(fig2B_out, device = "jpeg", unit = "px", width = 2160, height = 2160, plot = fig2B)
  
  vebcat("Successfully visualized highest spread", color = "funSuccess")
}

visualize_suitability <- function(rast, region, region.name, plot.show = F, verbose = F) {
  
  vebcat("Visualizing suiability plot", color = "funInit")
  
  if (!identical(crs(rast, proj = TRUE), crs(region, proj = TRUE))) {
    catn("Reprojecting to laea.")
    catn(crs(rast, proj = TRUE))
    catn(crs(region, proj = TRUE))
    catn(identical(crs(rast, proj = TRUE), crs(region, proj = TRUE)))
    rast <- project(rast, laea_crs, method = "bilinear")
  }
  
  catn("Acquiring min and max values.")
  prob_min <- 0.001
  prob_max <- where.max(rast[[1]])[[3]]
  region_ext <- ext(region)
  
  catn("Generating plot.")
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
  
  fig3_out <- "./outputs/visualize/plots/figure-3.jpeg"
  
  catn("Saving plot to:", colcat(fig3_out, color = "output"))
  
  ggsave(fig3_out, device = "jpeg", unit = "px", width = 2160, height = 2160, plot = fig3)
  
  vebcat("Suitability plot successfully visualized", color = "funSuccess")
}

visualize_richness <- function(dt, axis.x, axis.y, fill, group, plot.show = F, verbose = F) {
  vebcat("Visualizing composition plot", color = "funInit")
  
  # Remove NAs
  dt <- dt[!is.na(dt[[fill]]), ]
  
  # Reorder the bars
  dt[[axis.x]] <- as.factor(dt[[axis.x]])
  dt[[fill]] <- as.factor(dt[[fill]])
  
  #vebcat("Before sorting:", veb = verbose)
  #vebprint(levels(dt[[fill]]), veb = verbose)
  
  levels(dt[[fill]]) <- order_by_apg(levels(dt[[fill]]), by = fill, verbose = verbose)
  
  # Order the entire data table by the group parameter
  #dt <- dt[order(dt[[group]]), ]
  # Order the levels of the x axis by the sorted dt 
  #dt[[axis.x]] <- factor(dt[[axis.x]], levels = unique(dt[[axis.x]]))
  
  #vebcat("After sorting:", veb = verbose)
  #vebprint(levels(dt[[fill]]), veb = verbose)
  
  fig4 <- ggplot(dt, aes(x = dt[[axis.x]], y = dt[[axis.y]], fill = dt[[fill]] )) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_whitebox_d(palette = "muted") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      x = "Floristic Region", 
      y = "Relative Richness", 
      fill = paste0(toupper(substr(fill, 1, 1)), substr(fill, 2, nchar(fill)))
    )
  
  if (plot.show) print(fig4)
  
  fig4_out <- "./outputs/visualize/plots/figure-4.jpeg"
  
  catn("Saving plot to:", colcat(fig4_out, color = "output"))
  
  ggsave(fig4_out, plot = fig4, width = 3000, height = 2160, device = "jpeg", unit="px")
  
  vebcat("Composition plot successfully visualized", color = "funSuccess")
}

visualize_sankey <- function(dt, taxon, level, plot.show = F, verbose = F) {
  vebcat("Visualizing data in a sankey plot", color = "funInit")
  
  dt_sank <- copy(dt)
  
  #dt_sank <- dt_sank[!is.na(dt_sank[[taxon]])]
  #dt_sank <- dt_sank[complete.cases(dt_sank[[taxon]])]
  
  catn("Creating sankey plot.")
  print(unique(dt_sank$country))
  
  fig5 <- ggplot(data = dt_sank, aes(axis1 = dt_sank[[level]], axis2 = country, y = relativeRichness)) +
    scale_x_discrete(limits = c("Origin", "Destination"), expand = c(.1, .1)) +
    labs(
      y = paste0("Relative ", toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon)), " Richness"),
      #fill = paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon))),
      title = paste0("Relative ", toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon)), " Richness from origin region to Arctic region")
    ) +
    geom_flow() +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(vjust = 0.5, hjust = 0.5)
    )
  
  
  if (plot.show) print(fig5)
  
  fig5_out <- "./outputs/visualize/plots/figure-5.jpeg"
  
  catn("Saving plot to:", colcat(fig5_out, color = "output"))
  
  ggsave(fig5_out, device = "jpeg", unit = "px", width = 3840, height = 3500, plot = fig5)

  vebcat("Sankey plot successfully visualized", color = "funSuccess")
}

visualize_connections <- function(dt, taxon, level, plot.show = F, verbose = F) {
  vebcat("Visualizing data in a sankey plot", color = "funInit")
  
  wm <- get_world_map(projection = longlat_crs)
  
  dt_sank <- copy(dt)
  
  #dt_sank <- dt_sank[!is.na(dt_sank[[taxon]])]
  #dt_sank <- dt_sank[complete.cases(dt_sank[[taxon]])]
  
  catn("Creating sankey plot.")
  print(unique(dt_sank$country))
  
  fig5 <- ggplot(data = dt_sank, aes(axis1 = dt_sank[[level]], axis2 = country, y = relativeRichness)) +
    scale_x_discrete(limits = c("Origin", "Destination"), expand = c(.1, .1)) +
    labs(
      y = paste0("Relative ", toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon)), " Richness"),
      #fill = paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon))),
      title = paste0("Relative ", toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon)), " Richness from origin region to Arctic region")
    ) +
    geom_flow() +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(vjust = 0.5, hjust = 0.5)
    )
  
  
  if (plot.show) print(fig5)
  
  fig5_out <- "./outputs/visualize/plots/figure-5.jpeg"
  
  catn("Saving plot to:", colcat(fig5_out, color = "output"))
  
  ggsave(fig5_out, device = "jpeg", unit = "px", width = 3840, height = 3500, plot = fig5)
  
  vebcat("Sankey plot successfully visualized", color = "funSuccess")
}




