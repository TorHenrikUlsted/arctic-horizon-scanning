#########################
#     Frequency 1A      #
#########################

visualize_freqpoly <- function(spec.cells, region, region.name, vis.x, vis.color, vis.shade,vis.shade.name, vis.title = FALSE, vis.binwidth = 1, vis.x.scale = "log10", vis.peak.threshold = NULL, save.dir, save.device = "jpeg", save.unit = "px", plot.show = FALSE, verbose = FALSE) {
  
  create_dir_if(save.dir)
  
  # Figure 1A Whole CAVM
  vebcat("Creating histogram for the entire Region", color = "funInit")
  
  vebprint(class(spec.cells), verbose, "Input data class:")
  vebprint(head(spec.cells, 3), verbose, "Input data sample")
  vebprint(max(spec.cells[[vis.x]], na.rm = TRUE), verbose, paste0("Input max ", vis.x, ":"))
  
  vebcat("Creating plot 1A.", veb = verbose)
  
  vis_x <- spec.cells[[vis.x]]
  
  
  # Remove NA values
  spec.cells <- spec.cells[!is.na(vis_x), ]
  
  # Remove Inf values
  spec.cells <- spec.cells[vis_x != Inf, ]
  
  # Remove 0 values
  spec.cells <- spec.cells[vis_x != 0, ]
  
  vis_x <- spec.cells[[vis.x]]
  
  min_lim <- 1
  max_lim <- max(vis_x, na.rm = TRUE)
  vis_breaks <- c(min_lim, 5,10,50,100,500,1000, max_lim)
  
  if (!is.null(vis.x.scale)) {
    if (vis.x.scale == "log") {
      vis.binwidth <- (vis.binwidth/10)
    } else if (vis.x.scale == "log10") {
      vis.binwidth <- vis.binwidth/20
    } else if (vis.x.scale == "sqrt") {
      vis.binwidth <- vis.binwidth/4
      vis.peak.threshold <- 0.0002
    }
  } 
  catn("Using Binwidth:", vis.binwidth)
  
    #spec.cells, aes(x = vis_x)
  # After_Stat accesses ggplot processed data, count represents the count of data points in each bin of the histogram. Count is cell with that value. So we divide cell by total cells
  fig1A <- ggplot(spec.cells, aes(x = vis_x)) +
    geom_freqpoly(
      binwidth = vis.binwidth,
      aes(y = after_stat(count) / sum(after_stat(count)), color = after_stat(x)),
      na.rm = TRUE,
      show.legend = TRUE
    ) +
    scale_color_viridis_b(
      option = "B",
      guide = guide_legend(reverse = TRUE),
      breaks = vis_breaks
    ) +
    labs(
      x = paste0("Potential Species Richness ", if(!is.null(vis.x.scale)) {paste0("(", vis.x.scale, ")")} ), 
      y = paste0("Proportion of cells in ", region.name), 
      title = if (vis.title) paste0("Potential New Alien Species Richness in ", region.name),
      color = "Potential New Alien Richness"
      ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.01)) +
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
  
  if (!is.null(vis.x.scale)) {
   fig1A <- fig1A + scale_x_continuous(trans = vis.x.scale, limits = c(min_lim, max_lim), labels = function(x) format(x, big.mark = ",", digits = 1, scientific = FALSE))
  }
  
  if (plot.show) print(fig1A)
  
  # Calculate the peaks
  fig_data <- as.data.table(ggplot_build(fig1A)$data[[1]])

  peak_data <- find_peaks(fig_data, column = "y", threshold = vis.peak.threshold, verbose = verbose)
  
if (is.null(vis.x.scale)) {
  peak_data <- peak_data %>%
    mutate(
      cellCount = count,
      count = x
    )
} else {
  peak_data <- peak_data %>%
    mutate(
      cellCount = count,
      count = case_when(
        vis.x.scale == "log10" ~ 10^x,
        vis.x.scale == "sqrt" ~ x^2,
        vis.x.scale == "log" ~ exp(x),
        TRUE ~ x
      )
    )
}
  
  # Add text labels at the peaks
  fig1A_numbers <- fig1A + 
    geom_text_repel(
      data = peak_data, 
      aes(x = count, y = y, label = paste0("(", sub("0\\.", ".", round(y, 3)), ", ", round(count, 0), ")"))
  ) +
    geom_point(
      data = peak_data, 
      aes(x = count, y = y)
    )
  
  save_ggplot(
    save.plot = fig1A_numbers, 
    save.name = "figure-1A-descriptive",
    save.width = 3000, 
    save.height = 2160,
    save.dir = save.dir, 
    save.device = save.device,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  save_ggplot(
    save.plot = fig1A, 
    save.name = "figure-1A",
    save.width = 3000, 
    save.height = 2160,
    save.dir = save.dir, 
    save.device = save.device,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  vebcat("Frequency for the entire region successfully visualized", color = "funSuccess")
  
  
  #########################
  #     Frequency 1B      #
  #########################
  
  
  # Figure 1B different regions
  vebcat("Creating histogram for each region in the", region.name, color = "funInit")
  
  # Order by vis.color
  spec.cells <- spec.cells[order(spec.cells[[vis.color]]), ]
  spec.cells[[vis.color]] <- factor(spec.cells[[vis.color]], levels = unique(spec.cells[[vis.color]]))
  
  vis_x <- spec.cells[[vis.x]]
  vis_color <- spec.cells[[vis.color]]
  vis_shade <- spec.cells[[vis.shade]]
  
  # ADD shades of country for each floreg
  plcol <- as.factor(vis_color)
  plshd <- as.factor(vis_shade)
  vis_amax <- 1
  vis_amin <- 0.1
  
  vebprint(levels(plcol), verbose, text = "Color levels:")
  vebprint(levels(plshd), verbose, text = "Shade levels:")
  
  fig1B <- ggplot(spec.cells, aes(x = vis_x, color = plcol, alpha = plshd))  +
    geom_freqpoly(
      binwidth = 0.1, 
      aes(y =  after_stat(count) / sum(after_stat(count)))
      ) + 
    scale_color_viridis_d(guide = "legend", option = "B") +
    scale_alpha_discrete(guide = "none", range = c(vis_amin, vis_amax)) +
    labs(
      x = paste0("Potential Species Richness ", if(!is.null(vis.x.scale)) {paste0("(", vis.x.scale, ")")}), 
      y = paste0("Proportion of cells in the different regions of the ", region.name), 
      title = if (vis.title) paste0("Potential New Alien Species Richness in Different Regions of the ", region.name), 
      color = "Country", 
      alpha = vis.shade.name
    ) +
    scale_y_continuous(limits = c(0, NA), labels = function(x) format(x, big.mark = ",", scientific = FALSE) , breaks = seq(0, 1, by = 0.01)) +
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
  
  if (!is.null(vis.x.scale)) {
    fig1B <- fig1B + scale_x_continuous(trans = vis.x.scale, limits = c(min_lim, max_lim), labels = function(x) format(x, big.mark = ",", digits = 1, scientific = FALSE))
  }
  
  save_ggplot(
    save.plot = fig1B, 
    save.name = "figure-1B",
    save.width = 3000, 
    save.height = 2160,
    save.dir = save.dir, 
    save.device = save.device,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  vebcat("Frequency for each floristic region successfully visualized", color = "funSuccess")
}

#########################
#       Hotspots 2      #
#########################

visualize_hotspots <- function(rast, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.device = "jpeg", save.unit = "px", plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Potential Species hotspots", color = "funInit")
  
  catn("Checking crs.")
  rast <- check_crs(rast, projection, "near")
  region <- check_crs(region, projection, "near")

  catn("Getting min and max values.")
  min_lim <- 0
  max_lim <- where.max(rast)
  max_lim <- max(max_lim[, 3])
  vis_breaks <- c(1,5,10,50,100,200,400,700,900, max_lim)
  catn("Plotting hotspots.")
  
  fig2A <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = rast) +
    ggplot.filler(
      gradient = vis.gradient,
      scale.variable = "b",
      limits = c(min_lim, max_lim), 
      breaks = vis_breaks,
      guide = guide_legend(reverse = TRUE, title.position = "top", label.position = "bottom", nrow = 1),
      na.value = "transparent"
    ) + 
    labs(
      title = if (vis.title) paste0("Potential New Alien Species Hotspots in the ", region.name), 
      fill = "Potential Species Richness") +
    coord_sf(xlim = c(extent$xmin, extent$xmax), ylim = c(extent$ymin, extent$ymax)) +
    theme_minimal() +
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 12, face = "bold.italic"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom"
    )
  
  save_ggplot(
    save.plot = fig2A, 
    save.name = "figure-2A",
    save.width = 3050, 
    save.height = 3000,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  vebcat("Successfully visualized Potential Species hotspots", color = "funSuccess")
}

#########################
#    Distribution 3A    #
#########################

visualize_distributions <- function(rast, region, region.name, extent, projection, projection.method, vis.gradient ="viridis-b", vis.wrap = FALSE, vis.title = FALSE, save.dir, save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing species with highest distributions", color = "funInit")
  
  catn("Checking crs.")
  rast <- check_crs(rast, projection, projection.method, verbose = verbose)
  region <- check_crs(region, projection, projection.method, verbose = verbose)  
  
  catn("Plotting hotspots.")
  
  fig3A <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = rast) +
    ggplot.filler(
      gradient = vis.gradient,
      breaks = c(0, 1), 
      labels =  c("0", "1"),
      end = 0.6,
      guide = guide_legend(reverse = TRUE, title.position = "top", label.position = "bottom", nrow = 1),
      na.value = "transparent"
    ) + 
    labs(
      title = if (vis.title) paste0("New Aliens with Highest Potential Area of Occupancy in ", region.name), 
      fill = "Potential Area of occupancy") +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax), 
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 12, face = "bold.italic"),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom",
      strip.text = element_text(size = 10, face = "italic")
    )
  
  if (!is.null(vis.wrap)) {
    fig3A <- fig3A + facet_wrap(~lyr, nrow = vis.wrap, ncol = 3, labeller = label_wrap_gen(width = 20))
  }
  
  if (plot.save) {
    save_ggplot(
      save.plot = fig3A, 
      save.name = "figure-3A",
      save.width = 2700, 
      save.height = 3000,
      save.dir = save.dir, 
      save.device = save.device,
      save.unit = save.unit,
      vis.title = vis.title, 
      plot.show = plot.show, 
      verbose = verbose
    )
  }
  
  vebcat("Successfully visualized highest spread", color = "funSuccess")
  
  return(fig3A)
}

#########################
#    Suitability 3B     #
#########################

visualize_suitability <- function(stack, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.wrap = TRUE, vis.title = FALSE, save.dir, save.name = "figure-3B", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  
  vebcat("Visualizing suitability plot", color = "funInit")
  
  stack <- check_crs(stack, projection, "bilinear", verbose = verbose)
  region <- check_crs(region, projection, "bilinear", verbose = verbose)
  
  if (!identical(crs(stack, proj = TRUE), crs(region, proj = TRUE))) {
    catn("Reprojecting to laea.")
    catn(crs(stack, proj = TRUE))
    catn(crs(region, proj = TRUE))
    catn(identical(crs(stack, proj = TRUE), crs(region, proj = TRUE)))
    stack <- project(stack, laea_crs, method = "bilinear")
  }
  
  catn("Acquiring min and max values.")
  min_lim <- 0
  max_lim <- where.max(stack[[1]])[1,][[3]] # Get the first row and last item "value"
  region_ext <- ext(region)
  vis_breaks <- seq(min_lim, max_lim, length.out = 50)
  
  catn("Generating plot.")
  fig3B <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = stack) +
    ggplot.filler(
      gradient = vis.gradient,
      limits = c(min_lim, max_lim), 
      breaks = c(seq(min_lim, max_lim, by = 50), max_lim), 
      guide = guide_legend(reverse = TRUE, title.position = "top", label.position = "bottom", nrow = 1),
      na.value = "transparent"
    ) + 
    labs(
      title = if (vis.title) paste0("Potential New Alien Climatic suitability in ", region.name), 
      fill = paste0("Potential Climatic Suitability", if (!is.null(vis.unit)) {paste0(" (", vis.unit, ")")} ), 
      color = "Country"
    ) +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax), 
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 12, face = "bold.italic"),
      axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom",
      strip.text = element_text(size = 10, face = "italic")
      )
  
  if (!is.null(vis.wrap)) {
    fig3B <- fig3B + facet_wrap(~lyr, nrow = vis.wrap, ncol = 3, labeller = label_wrap_gen(width = 20))
  } 
  
  if (plot.show) print(fig3B)
  
  if (plot.save) {
    save_ggplot(
      save.plot = fig3B, 
      save.name = save.name,
      save.width = 2700, 
      save.height = 3000,
      save.dir = save.dir,
      save.device = save.device,
      save.unit = save.unit,
      vis.title = vis.title,
      plot.show = plot.show, 
      verbose = verbose
    )
  }
  
  
  vebcat("Suitability plot successfully visualized", color = "funSuccess")
  
  return(fig3B)
}

##########################################
#         suitability Units 3C           #
##########################################

visualize_suit_units <- function(stack.mean, stack.median, stack.max, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.wrap = NULL, vis.title = FALSE, save.dir, save.name = "figure-3C", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  
  n <- vis.wrap * 3
  
  stack_mean <- stack.mean[[1:n]]
  stack_median <- stack.median[[1:n]]
  stack_max <- stack.max[[1:n]]
  
  fig_mean <- visualize_suitability(
    stack = stack_mean,
    region = world_map, 
    region.name = region.name,
    extent = region_ext,
    projection = laea_crs,
    vis.gradient = vis.gradient,
    vis.wrap = vis.wrap,
    save.dir = plot_dir,
    save.device = save.device,
    plot.save = FALSE,
    verbose = verbose
  ) 
  
  fig_mean + theme(legend.position = "none")
  
  fig_median <- visualize_suitability(
    stack = stack_median,
    region = world_map, 
    region.name = region.name,
    extent = region_ext,
    projection = laea_crs,
    vis.gradient = vis.gradient,
    vis.wrap = vis.wrap,
    save.dir = plot_dir,
    save.device = save.device,
    plot.save = FALSE,
    verbose = verbose
  )
  
  fig_median + theme(legend.position = "none")
  
  fig_max <- visualize_suitability(
    stack = stack_max,
    region = world_map, 
    region.name = region.name,
    extent = region_ext,
    projection = laea_crs,
    vis.gradient = vis.gradient,
    vis.wrap = vis.wrap,
    save.dir = plot_dir,
    save.device = save.device,
    plot.save = FALSE,
    verbose = verbose
  )
  
  fig3C <- ggarrange(
    fig_mean, fig_median, fig_max,
    labels = c("Suitability (mean)\n\n\n", "Suitability (median)\n\n\n", "\nSuitability (max)\n\n\n"),
    ncol = 1,
    nrow = n
  )
  
  save_ggplot(
    save.plot = fig3C, 
    save.name = "figure-3C",
    save.width = 3200, 
    save.height = 4000,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
}

##########################################
#     Distribution + suitability 3D      #
##########################################

visualize_dist_suit <- function(stack.distribution, stack.suitability, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.wrap = NULL, vis.title = FALSE, save.dir, save.name = "figure-3D", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  
  n <- vis.wrap * 3
  
  dist_raster <- stack.distribution[[1:n]]
  suit_raster <- stack.suitability[[1:n]]
  
    fig_dist <- visualize_distributions(
      rast = dist_raster,
      region = world_map,
      region.name = region.name,
      extent = region_ext,
      projection  = out_projection,
      vis.gradient = vis.gradient,
      vis.wrap = vis.wrap,
      save.dir = plot_dir,
      plot.save = FALSE,
      verbose = verbose
    )
    
    fig_suit <- visualize_suitability(
      stack = suit_raster,
      region = world_map, 
      region.name = region.name,
      extent = region_ext,
      projection = laea_crs,
      vis.gradient = vis.gradient,
      vis.wrap = vis.wrap,
      save.dir = plot_dir,
      save.device = save.device,
      plot.save = FALSE,
      verbose = verbose
    )
    
  fig3D <- ggarrange(
    fig_dist, fig_suit,
    labels = c("Distribution", "Suitability"),
    ncol = 1,
    nrow = 2
  )
  
  save_ggplot(
    save.plot = fig3D, 
    save.name = "figure-3D",
    save.width = 2700, 
    save.height = 3000,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
}


visualize_richness <- function(dt, region.name, vis.x, vis.x.sort, vis.y, vis.fill, vis.group, vis.gradient, vis.title = FALSE, save.dir, save.device = "jpeg", save.unit = "px", plot.show = F, verbose = F) {
  vebcat("Visualizing composition plot", color = "funInit")
  
  # Remove NAs
  #dt <- dt[!is.na(dt[[vis.fill]]), ]
  
  print(dt)
  
  vebprint(length(unique(dt[[vis.x]])), text = "Length of Unique Sub Regions:")
  # 
  # # Order dt by vis.x.sort column
  # dt <- dt[order(dt[[vis.x.sort]])]
  
  print(dt)
  
  lvl_order <- unique(dt[[vis.x]])
  
  dt[[vis.x]] <- factor(dt[[vis.x]], levels = lvl_order)
  
  # Reorder the bars
  dt[[vis.fill]] <- factor(dt[[vis.fill]])
  
  levels(dt[[vis.fill]]) <- order_by_apg(levels(dt[[vis.fill]]), by = vis.fill, verbose = verbose)
  
  print(levels(dt[[vis.fill]]))
  
  fig4 <- ggplot(dt, aes(x = dt[[vis.x]], y = dt[[vis.y]], fill = dt[[vis.fill]] )) +
    geom_bar(stat = "identity", position = "stack") +
    ggplot.filler(
      gradient = vis.gradient,
      scale.variable = "d",
      guide = guide_legend(reverse = FALSE),
      na.value = "transparent"
    ) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(
      x = "Floristic Region",
      y = "Potential Relative Species Richness",
      title = paste0("Potential New Alien ", paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill))), " Composition in ", region.name),
      fill = paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill)))
    ) +
    theme_minimal() + 
    theme(
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 12, face = "bold.italic"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10, hjust = 0.5),
    )
  
  if (plot.show) print(fig4)
  
  save_ggplot(
    save.plot = fig4, 
    save.name = "figure-4",
    save.width = 5000, 
    save.height = 4000,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
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

visualize_connections <- function(dt, taxon, region.name, vis.gradient, vis.title = FALSE, save.dir, save.name = "figure-3D", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing data in a sankey plot", color = "funInit")
  
  wm <- get_world_map(projection = mollweide_crs)
  origin_points <- get_con_points(dt, "mollweide", "originMeanLong", "originMeanLat", verbose = verbose)
  dest_points <- get_con_points(dt, "mollweide", "subRegionLong", "subRegionLat", verbose = verbose)
  
  origin_points <- origin_points[!is.na(origin_points)]
  dest_points <- dest_points[!is.na(dest_points)]
  
  origin <- terra::as.data.frame(origin_points, geom = "xy")
  origin <- as.data.table(origin)
  setnames(origin, "x", "originX")
  setnames(origin, "y", "originY")
  
  dest <- terra::as.data.frame(dest_points, geom = "xy")
  dest <- as.data.table(dest)
  setnames(dest, "x", "destX")
  setnames(dest, "y", "destY")
  
  merged_dt <- merge(dest, origin, by = "cleanName")
  
  print(merged_dt)
           
  catn("Creating connections plot.")
  
  fig5 <- ggplot(data = merged_dt) +
    geom_spatvector(data = wm) +
    # geom_spatvector(data = origin_points, aes(color = "Origin")) +
    # geom_spatvector(data = dest_points, aes(color = "Destination")) +
    geom_point(aes(x = originX, y = originY, color = "Origin")) +
    geom_point(aes(x = destX, y = destY, color = "Destination")) +
    geom_segment(aes(x = originX, y = originY, xend = destX, yend = destY)) +
    # coord_sf(crs = mollweide_crs)+
    labs(
      x = "Longitude",
      y = "Latitude",
      color = c("Regions"),
      color = paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon))),
      title = paste0("Potential New Alien ", paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon))), " Richness from origin Country to Arctic Florsitic Province")
    ) +
    scale_size_continuous(range = c(1, 1.5)) +
    labs(color = "Order Count") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(vjust = 0.5, hjust = 0.5)
    )
    
  if (plot.show) print(fig5)
  
  save_ggplot(
    save.plot = fig5, 
    save.name = "figure-5",
    save.width = 4000, 
    save.height = 3000,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  vebcat("Connection plot successfully visualized", color = "funSuccess")
}




