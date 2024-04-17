visualize_freqpoly <- function(sp_cells, region, region.name, vis.x, vis.color, vis.shade,vis.shade.name, vis.title = FALSE, plot.show = FALSE, verbose = FALSE) {
  
  create_dir_if("./outputs/visualize/plots")
  
  # Figure 1A Whole CAVM
  vebcat("Creating histogram for the entire Region", color = "funInit")
  
  print(class(sp_cells))
  print(head(sp_cells, 3))
  
  vebcat("Creating plot 1A.", veb = verbose)
  
  vis_x <- sp_cells[[vis.x]]
  #sp_cells[, richness := ifelse(richness == 0, NA, richness)]
  # Remove NA values
  sp_cells <- sp_cells[!is.na(vis_x), ]
  
  # Remove Inf values
  sp_cells <- sp_cells[vis_x != Inf, ]
  
  # Remove 0 values
  sp_cells <- sp_cells[vis_x != 0, ]
  
  vis_x <- sp_cells[[vis.x]]
  
  min_lim <- 0
  max_lim <- max(vis_x)
  vis_breaks <- c(1,5,10,50,100,500,1000, max_lim)
  
    #sp_cells, aes(x = vis_x)
  # After_Stat accesses ggplot processed data, count represents the count of data points in each bin of the histogram. Count is cell with that value. So we divide cell by total cells
  fig1A <- ggplot(sp_cells, aes(x = vis_x)) +
    geom_freqpoly(
      binwidth = 0.1, 
      aes(y = after_stat(count) / sum(after_stat(count)), color = after_stat(x)),
      na.rm = TRUE,
      show.legend = TRUE
    ) +
    scale_color_viridis_b(
      option = "B",
      guide = guide_legend(reverse = TRUE),
      limits = (c(min_lim, max_lim)),
      breaks = vis_breaks
    ) +
    labs(
      x = "Potential Species Richness (log10)", 
      y = paste0("Proportion of cells in the ", region.name), 
      title = if (vis.title) paste0("Potential New Alien Species Richness in the ", region.name),
      color = "Potential Species Richness"
      ) +
    scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme_minimal() + 
    theme(
      plot.title = element_text(
        color = "black", 
        vjust = -0.5, 
        hjust = 0.5, 
        size = 14, 
        face = "bold.italic"
      ),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(
        hjust = 0.5
      )
    ) 
  
  if (plot.show) print(fig1A)
  
  # Calculate the peaks
  fig_data <- as.data.table(ggplot_build(fig1A)$data[[1]])

  peak_data <- find_peaks(fig_data, column = "y", threshold = NULL, verbose = verbose)
    
  peak_data <- peak_data %>% 
    mutate(cellCount = count, count = round(10^x))
  
  # Add text labels at the peaks
  fig1A <- fig1A + geom_text(
    data = peak_data, 
    aes(x = count, y = y, label = count),
    vjust = -0.5,
    hjust = 0.5
  )
  
  ggsave("./outputs/visualize/plots/figure-1A.jpeg", device = "jpeg", unit = "px", width = 3000, height = 2160, fig1A)
  
  vebcat("Frequency for the entire region successfully visualized", color = "funSuccess")
  
  
# ----------- Sub regions -------------- #
  
  
  # Figure 1B different regions
  vebcat("Creating histogram for each region in the", region.name, color = "funInit")
  
  # Order by vis.color
  sp_cells <- sp_cells[order(sp_cells[[vis.color]]), ]
  sp_cells[[vis.color]] <- factor(sp_cells[[vis.color]], levels = unique(sp_cells[[vis.color]]))
  
  vis_x <- sp_cells[[vis.x]]
  vis_color <- sp_cells[[vis.color]]
  vis_shade <- sp_cells[[vis.shade]]
  
  # ADD shades of country for each floreg
  plcol <- as.factor(vis_color)
  plshd <- as.factor(vis_shade)
  vis_amax <- 1
  vis_amin <- 0.1
  
  vebprint(levels(plcol), verbose, text = "Color levels:")
  vebprint(levels(plshd), verbose, text = "Shade levels:")
  
  fig1B <- ggplot(sp_cells, aes(x = vis_x, color = plcol, alpha = plshd))  +
    geom_freqpoly(
      binwidth = 0.1, 
      aes(y =  after_stat(count) / sum(after_stat(count)))
      ) + 
    scale_color_viridis_d(guide = "legend", option = "B") +
    scale_alpha_discrete(guide = "none", range = c(vis_amin, vis_amax)) +
    labs(
      x = "Potential Species Richness (log10)", 
      y = paste0("Proportion of cells in the different regions of the ", region.name), 
      title = if (vis.title) paste0("Potential New Alien Species Richness in Different Regions of the ", region.name), 
      color = "Country", 
      alpha = vis.shade.name
    ) +
    scale_y_continuous(limits = c(0, NA), labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    scale_x_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    theme_minimal() + 
    theme(
      plot.title = element_text(
        color = "black", 
        vjust = -0.5, 
        hjust = 0,
        size = 14, 
        face = "bold.italic"
      )
    )
  
  if (plot.show) print(fig1B)
  
  ggsave("./outputs/visualize/plots/figure-1B.jpeg",device = "jpeg", unit = "px", width = 3000, height = 2160, fig1B)
  
  vebcat("Frequency for each floristic region successfully visualized", color = "funSuccess")
}

# ---------- Hot spots ----------- #

visualize_hotspots <- function(rast, region, region.name, extent, projection, projection.method, vis.title = FALSE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Potential Species hotspots", color = "funInit")
  
  catn("Checking crs.")
  rast <- check_crs(rast, projection, projection.method)
  region <- check_crs(region, projection, "near")
  
  catn("Removing 0 values.")
  # Remove 0 values
  #rast_0 <- ifel(rast == 0, rast, NA)
  #rast <- ifel(rast == 0, NA, rast)

  catn("Getting min and max values.")
  min_lim <- 0
  max_lim <- where.max(rast)
  max_lim <- max(max_lim[, 3])
  vis_breaks <- c(1,5,10,50,100,500,1000, max_lim)
  catn("Plotting hotspots.")
  
  fig2A <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = rast) +
    scale_fill_viridis_b(
      option = "B", 
      guide = guide_legend(reverse = TRUE),
      limits = (c(min_lim, max_lim)),
      breaks = vis_breaks,
      labels = function(x) format((x), big.mark = ",", scientific = FALSE, digits = 2),
      na.value = "transparent"
    ) +
    labs(title = if (vis.title) paste0("Potential New Alien Species Hotspots in the ", region.name), fill = "Potential Species Richness") +
    coord_sf(xlim = c(extent$xmin, extent$xmax), ylim = c(extent$ymin, extent$ymax)) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 12, face = "bold.italic"),
      legend.title = element_text(size = 8, hjust = 0),
      legend.text = element_text(size = 8)
    )
  
  if (plot.show) print(fig2A)
  
  catn("Saving plot.")
  ggsave("./outputs/visualize/plots/figure-2A.jpeg", device = "jpeg", unit = "px", width = 2700, height = 2160, plot = fig2A)
  
  vebcat("Successfully visualized Potential Species hotspots", color = "funSuccess")
}

# ---------- Hot spots 2B ----------- #

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
    scale_fill_whitebox_c(
      palette = "muted", 
      breaks = c(0, 1), 
      labels = c("No Climatic Overlap", "Climatic Overlap"), 
      guide = guide_legend(reverse = TRUE)
    ) +
    labs(
      title = if (vis.title) paste0("Species with Highest Potential Spread in the", region.name), 
      fill = "Overlap Value") +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax), 
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, size = 12, face = "bold.italic"),
      legend.title = element_text(size = 5),
      legend.text = element_text(size = 8)
    )
  
  if (plot.show) print(fig2B)
  
  fig2B_out <- "./outputs/visualize/plots/figure-2B.jpeg"
  
  catn("Saving plot to:", colcat(fig2B_out, color = "output"))
  ggsave(fig2B_out, device = "jpeg", unit = "px", width = 2700, height = 2160, plot = fig2B)
  
  vebcat("Successfully visualized highest spread", color = "funSuccess")
}

# ---------- Suitability ----------- #

visualize_suitability <- function(rast, region, region.name, vis.title = FALSE, plot.show = F, verbose = F) {
  
  vebcat("Visualizing suiability plot", color = "funInit")
  
  rast <- check_crs(rast, projection, projection.method)
  region <- check_crs(region, projection, projection.method)
  
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
    labs(
      title = if (vis.title) paste0("Predicted species distribution in ", region.name), 
      fill = "Potential Species Suitability", 
      color = "Country"
    ) +
    coord_sf(xlim = c(region_ext$xmin, region_ext$xmax), ylim = c(region_ext$ymin, region_ext$ymax)) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 10, face = "bold.italic"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      legend.text = element_text(size = 5),
      strip.text = element_text(size = 10)
      )
  
  if (plot.show) print(fig3)
  
  fig3_out <- "./outputs/visualize/plots/figure-3.jpeg"
  
  catn("Saving plot to:", colcat(fig3_out, color = "output"))
  
  ggsave(fig3_out, device = "jpeg", unit = "px", width = 2160, height = 2160, plot = fig3)
  
  vebcat("Suitability plot successfully visualized", color = "funSuccess")
}

visualize_richness <- function(dt, vis.x, vis.x.sort, vis.y, vis.fill, vis.group, plot.show = F, verbose = F) {
  vebcat("Visualizing composition plot", color = "funInit")
  
  # Remove NAs
  dt <- dt[!is.na(dt[[vis.fill]]), ]
  
  vebprint(length(unique(dt[[vis.x]])), verbose, "Length of Unique Sub Regions:")
  
  # Order dt by vis.x.sort column
  dt <- dt[order(dt[[vis.x.sort]])]
  
  lvl_order <- unique(dt[[vis.x]])
  
  dt[[vis.x]] <- factor(dt[[vis.x]], levels = lvl_order)
  
  # Reorder the bars
  dt[[vis.fill]] <- factor(dt[[vis.fill]])
  
  levels(dt[[vis.fill]]) <- order_by_apg(levels(dt[[vis.fill]]), by = vis.fill, verbose = verbose)
  
  fig4 <- ggplot(dt, aes(x = dt[[vis.x]], y = dt[[vis.y]], fill = dt[[vis.fill]] )) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_whitebox_d(palette = "muted") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      x = "Floristic Region", 
      y = "Potential Relative Richness", 
      fill = paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill)))
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

visualize_connections <- function(dt, plot.show = F, verbose = F) {
  vebcat("Visualizing data in a sankey plot", color = "funInit")
  
  wm <- get_world_map(projection = longlat_crs)
  
  con_dt <- copy(dt)
  
  # Handle endcoords and floristicprovince if needed
  
  con_dt <- con_dt[!is.na(con_dt$country), ]
  
  # then calc mean_lat and _mean_long for each countryCode
  # Count the number of unique species for each countryCode
  dt_means <- con_dt[, .(meanLongOrig = mean(meanLong, na.rm = TRUE), 
                         meanLatOrig = mean(meanLat, na.rm = TRUE),
                         meanLongDest = mean(longDest, na.rm = TRUE),
                         meanLatDest = mean(latDest, na.rm = TRUE),
                     species_count = uniqueN(cleanName)), 
                 by = .("country", "floristicProvince")]
  # Get endLong and endLat
  
  catn("Creating connections plot.")

  ggplot(data = dt_means) +
    geom_spatvector(data = world_map) +
    scale_x_continuous(name = "Longitude", breaks = seq(-180, 180, 30)) +
    scale_y_continuous(name = "Latitude", breaks = seq(-90, 90, 15)) +
    geom_point(aes(x = meanLong, y = meanLat, color = "Origin")) +
    geom_point(aes(x = endLong, y = endLat, color = "Destination")) +
    labs(
      color = c("Regions"),
      #    color = paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon)), " Count"),
      #fill = paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon))),
      title = "Richness from origin region to Arctic region"
    ) +
    new_scale_color() +
    geom_segment(aes(x = meanLong, y = meanLat, xend = endLong, yend = endLat, color = species_count)) +
    scale_size_continuous(range = c(1, 1.5)) +
    scale_colour_viridis_c(
      option = "B", 
      guide = guide_legend(reverse = TRUE),
      #limits = (c(min_lim, max_lim)),
      #breaks = vis_breaks,
      labels = function(x) format((x), big.mark = ",", scientific = FALSE, digits = 2),
      na.value = "transparent"
    ) +
    labs(color = "Order Count") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(vjust = 0.5, hjust = 0.5)
    )
    
  if (plot.show) print(fig6)
  
  fig6_out <- "./outputs/visualize/plots/figure-6.jpeg"
  
  catn("Saving plot to:", colcat(fig6_out, color = "output"))
  
  ggsave(fig6_out, device = "jpeg", unit = "px", width = 3840, height = 3500, plot = fig6)
  
  vebcat("Connection plot successfully visualized", color = "funSuccess")
}




