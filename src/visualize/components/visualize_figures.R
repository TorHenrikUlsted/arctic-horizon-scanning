#########################
#     Frequency 1A      #
#########################

visualize_freqpoly <- function(spec.cells, region, region.name, vis.x, vis.color, vis.shade,vis.shade.name, vis.gradient = "viridis-B", vis.title = FALSE, vis.binwidth = 1, vis.x.scale = "log", vis.peak.threshold = NULL, save.dir, save.device = "jpeg", save.unit = "px", plot.show = FALSE, verbose = FALSE) {
  
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
  #spec.cells <- spec.cells[vis_x != 0, ]
  
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
    ggplot.filler(
      gradient = gradient.config$vis.gradient,
      scale.type = "color-b",
      guide = gradient.config$guide,
      breaks = vis_breaks,
      labels = function(x) sprintf("%.0f", round(x, 0)),
      na.value = gradient.config$na.value
    ) +
    labs(
      x = paste0("Potential Species Richness ", if(!is.null(vis.x.scale)) {paste0("(", vis.x.scale, ")")} ), 
      y = paste0("Proportion of cells in ", region.name), 
      title = if (vis.title) paste0("Potential New Alien Species Richness in ", region.name),
      color = "Potential Species Richness"
      ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.01)) +
    theme_minimal() + 
    theme.config +
    theme(
      axis.title.x = element_text(color = "#575757"),
      axis.title.y = element_text(color = "#575757"),
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
    ggplot.filler(
      gradient = gradient.config$vis.gradient,
      guide = gradient.config$guide,
      scale.type = "color-d",
      breaks = vis_breaks,
      labels = function(x) sprintf("%.0f", round(x, 0)),
      na.value = gradient.config$na.value
    ) +
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
    theme.config +
    theme(
      axis.title.x = element_text(color = "#575757"),
      axis.title.y = element_text(color = "#575757"),
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

visualize_hotspots <- function(raster, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.device = "jpeg", save.unit = "px", plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Potential Species hotspots", color = "funInit")
  # 
  # catn("Checking crs.")
  # #raster <- check_crs(raster, projection, "near")
  # region <- check_crs(region, projection, "near")

  catn("Getting min and max values.")
  min_lim <- 0
  max_lim <- where.max(raster)
  max_lim <- max(max_lim[, 3])
  vis_breaks <- c(1,5,10,50,100,200,400,700,900, max_lim)
  catn("Plotting hotspots.")
  
  fig2A <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = raster) +
    ggplot.filler(
      gradient = gradient.config$vis.gradient,
      guide = gradient.config$guide,
      scale.type = "fill-b",
      limits = c(min_lim, max_lim), 
      breaks = vis_breaks,
      labels = function(x) sprintf("%.0f", round(x, 0)),
      na.value = gradient.config$na.value
    ) + 
    labs(
      title = if (vis.title) paste0("Potential New Alien Species Hotspots in the ", region.name), 
      fill = "Potential Species Richness") +
    coord_sf(xlim = c(extent$xmin, extent$xmax), ylim = c(extent$ymin, extent$ymax)) +
    theme_minimal() +
    theme.config
  
  save_ggplot(
    save.plot = fig2A, 
    save.name = "figure-2",
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

visualize_paoo <- function(rast, region, region.name, extent, projection, projection.method, vis.gradient ="viridis-b", vis.wrap = FALSE, vis.title = FALSE, save.dir, save.device = "jpeg", save.unit = "px", return = FALSE, plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing species with highest Potential Area of occupancy", color = "funInit")
  
  catn("Checking crs.")
  rast <- check_crs(rast, projection, projection.method, verbose = verbose)
  region <- check_crs(region, projection, projection.method, verbose = verbose)  
  
  catn("Plotting hotspots.")
  
  fig3A <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = rast) +
    ggplot.filler(
      gradient = gradient.config$vis.gradient,
      guide = gradient.config$guide,
      scale.type = "fill-c",
      breaks = c(0, 1), 
      labels =  c("0", "1"),
      end = 0.6,
      na.value = gradient.config$na.value
    ) + 
    labs(
      title = if (vis.title) paste0("Potential New Aliens with Highest Potential Area of Occupancy in ", region.name), 
      fill = "Potential Area of occupancy") +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax), 
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() + 
    theme.config +
    theme(
      strip.text = element_text(size = 16, face = "italic")
    )
  
  if (!is.null(vis.wrap)) {
    fig3A <- fig3A + facet_wrap(~lyr, nrow = vis.wrap, ncol = 3, labeller = label_wrap_gen(width = 50))
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
  
  vebcat("Successfully visualized highest Potential Area of Occupancy", color = "funSuccess")
  
  if(return) return(fig3A)
}

#########################
#    Suitability 3B     #
#########################

visualize_suitability <- function(stack, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.wrap = TRUE, vis.title = FALSE, save.dir, save.name = "figure-3B", save.device = "jpeg", save.unit = "px", plot.save = TRUE, return = FALSE, plot.show = FALSE, verbose = FALSE) {
  
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
  min_lim <- 0.00
  max_lim <- where.max(stack[[1]])[1,][[3]] # Get the first row and last item "value"
  region_ext <- ext(region)
  vis_breaks <- seq(min_lim, max_lim, length.out = 8)
  
  catn("Generating plot.")
  fig3B <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = stack) +
    ggplot.filler(
      gradient = gradient.config$vis.gradient,
      guide = gradient.config$guide,
      scale.type = "fill-c",
      limits = c(min_lim, max_lim), 
      breaks = vis_breaks,
      labels = function(x) sprintf("%.2f", round(x, 2)),
      na.value = gradient.config$na.value
    ) + 
    labs(
      title = if (vis.title) paste0("Potential New Alien Climatic suitability in ", region.name), 
      fill = paste0("Potential Climatic Suitability"),# if (!is.null(vis.unit)) {paste0(" (", vis.unit, ")")} ), 
      color = "Country"
    ) +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax), 
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() + 
    theme.config +
    theme(
      strip.text = element_text(size = 16, face = "italic")
      )
  
  if (!is.null(vis.wrap)) {
    fig3B <- fig3B + facet_wrap(~lyr, nrow = vis.wrap, ncol = 3, labeller = label_wrap_gen(width = 50))
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
  
  if(return) return(fig3B)
}

##########################################
#         suitability Units 3C           #
##########################################

visualize_suit_units <- function(stack.mean, stack.median, stack.max, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.wrap = NULL, vis.title = FALSE, save.dir, save.name = "figure-3C", save.device = "jpeg", save.unit = "px", return = TRUE, plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  
  n <- vis.wrap * 3
  
  stack_mean <- stack.mean[[1:n]]
  stack_median <- stack.median[[1:n]]
  stack_max <- stack.max[[1:n]]
  
  fig_mean <- visualize_suitability(
    stack = stack_mean,
    region = region, 
    region.name = region.name,
    extent = extent,
    projection = projection,
    vis.wrap = vis.wrap,
    save.dir = save.dir,
    save.device = save.device,
    return = return,
    plot.save = FALSE,
    plot.show = plot.show,
    verbose = verbose
  ) 
  
  fig_mean <- fig_mean + theme(legend.position = "none") + ggtitle("Potential Mean Suitability")
  
  fig_median <- visualize_suitability(
    stack = stack_median,
    region = region, 
    region.name = region.name,
    extent = extent,
    projection = projection,
    vis.wrap = vis.wrap,
    save.dir = save.dir,
    save.device = save.device,
    return = return,
    plot.save = FALSE,
    plot.show = plot.show,
    verbose = verbose
  )
  
  fig_median <- fig_median + theme(legend.position = "none") + ggtitle("Potential Median Suitability")
  
  fig_max <- visualize_suitability(
    stack = stack_max,
    region = region, 
    region.name = region.name,
    extent = extent,
    projection = projection,
    vis.wrap = vis.wrap,
    save.dir = save.dir,
    save.device = save.device,
    return = return,
    plot.save = FALSE,
    verbose = verbose
  )
  
  fig_max <- fig_max + ggtitle("Potential Max Suitability")
  
  catn("Arranging plots.")
  
  fig3C <- grid.arrange(
    fig_mean, fig_median, fig_max,
    ncol = 1,  # One column for the first row
    heights = c(1, 1, 1)  # Equal heights for all rows
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

visualize_dist_suit <- function(stack.paoo, stack.suitability, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.wrap = NULL, vis.title = FALSE, save.dir, save.name = "figure-3D", save.device = "jpeg", save.unit = "px", return = TRUE, plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  
  n <- vis.wrap * 3
  
  paoo_raster <- stack.paoo[[1:n]]
  suit_raster <- stack.suitability[[1:n]]
  
    fig_paoo <- visualize_paoo(
      rast = paoo_raster,
      region = region,
      region.name = region.name,
      extent = extent,
      projection  = projection,
      vis.wrap = vis.wrap,
      save.dir = save.dir,
      return = return,
      plot.save = FALSE,
      verbose = verbose
    )
    
    fig_paoo <- fig_paoo + ggtitle("Potential Area of Occupancy")
    
    fig_suit <- visualize_suitability(
      stack = suit_raster,
      region = region, 
      region.name = region.name,
      extent = extent,
      projection = projection,
      vis.wrap = vis.wrap,
      save.dir = save.dir,
      save.device = save.device,
      return = return,
      plot.save = FALSE,
      verbose = verbose
    )
    
    fig_suit <- fig_suit + ggtitle("Potential Suitability")
    
    catn("Arranging plots.")
    
    fig3D <- grid.arrange(
      fig_paoo, fig_suit,
      ncol = 1,  # One column for the first row
      nrow = 2   # Two rows
    )
    
  save_ggplot(
    save.plot = fig3D, 
    save.name = save.name,
    save.width = 3200, 
    save.height = 3200,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title,
    plot.show = plot.show, 
    verbose = verbose
  )
}

##########################
#   Composition 4A & B   #
##########################

visualize_composition <- function(dt, region.name, vis.x, vis.x.sort, vis.y, vis.fill, vis.group, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.device = "jpeg", save.unit = "px", plot.show = F, verbose = F) {
  vebcat("Visualizing composition plot", color = "funInit")
  
  dt_copy <- copy(dt)
  
  vebcat("Number of Sub Regions:", highcat(length(unique(dt_copy[[vis.x]]))))
  # 
  # # Order dt by vis.x.sort column
  #setorder(dt, get(vis.x.sort))
  dt_copy <- dt_copy[order(dt_copy[[vis.x.sort]])]
  
  lvl_order <- unique(dt_copy[[vis.x]])
  
  dt_copy[[vis.x]] <- factor(dt_copy[[vis.x]], levels = lvl_order)
  
  # Reorder the bars
  dt_copy[[vis.fill]] <- factor(dt_copy[[vis.fill]])
  
  vebprint(levels(dt_copy[[vis.fill]]), verbose, "Order levels:")
  
  dt_copy <- get_order_group(dt_copy, verbose = verbose)
  
  vebprint(levels(dt_copy[[vis.fill]]), verbose, "Final Levels:")
  
  p1 <-  geom_bar(stat = "identity", position = "stack")
  
  # Save config
  saved_config_params <- gradient.config$guide$params
  gradient.config$guide$params$nrow <- NULL
  gradient.config$guide$params$ncol <- 1
  
  p2 <- ggplot.filler(
    gradient = gradient.config$vis.gradient,
    guide = guide_legend(ncol = 2),
    scale.type = "fill-d",
    na.value = gradient.config$na.value
  )
  p4 <- labs(
    x = "Floristic Province",
    y = paste0("Potential Relative ", 
               paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill))), 
               " Richness"
              ),
  title = if (vis.title) paste0(
    "Potential New Alien ", 
    paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill))), 
    " Composition in ", 
    region.name
  ),
    fill = paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill)))
  )
  p5 <- theme_minimal()
  p6 <- theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
    legend.position = "right"
  )
  
  fig4A <- ggplot(dt_copy, aes(x = get(vis.x), y = get(vis.y), fill = get(vis.fill))) + 
    p1+p2+p4+p5+theme.config+p6
  
  if (plot.show) print(fig4A)
  
  save_ggplot(
    save.plot = fig4A, 
    save.name = "figure-4A",
    save.width = 3700, 
    save.height = 3700,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  fig4B <- ggplot(dt_copy, aes(x = get(vis.x), y = get(vis.y), fill = get(vis.group))) + 
    p1 +
    ggplot.filler(
    gradient = vis.gradient,
    scale.type = "fill-d",
    guide = guide_legend(ncol = 1),
    na.value = gradient.config$na.value
  ) +
    p4+p5+theme.config+p6
  
  if (plot.show) print(fig4B)
  
  save_ggplot(
    save.plot = fig4B, 
    save.name = "figure-4B",
    save.width = 3700, 
    save.height = 3300,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  # reset config
  gradient.config$guide$params <- saved_config_params
  
  vebcat("Composition plot successfully visualized", color = "funSuccess")
}

#######################
#    Connections 5    #
#######################

visualize_connections <- function(dt, taxon, region.name, subregion.name, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-3D", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Connections map", color = "funInit")
  
  wm <- get_world_map(projection = mollweide_crs)
  sub_dt <- copy(dt)
  #sub_dt <- dt[, nLines := uniqueN(get(taxon), na.rm = TRUE), by = .(originCountry, subRegionName)]
  
  origin_subset <- sub_dt[, .(get(taxon), originMeanLong, originMeanLat, connections)]
  setnames(origin_subset, "V1", taxon)
  origin_subset <- unique(origin_subset, by = c(taxon, "originMeanLong", "originMeanLat"))
  origin_points <- get_con_points(origin_subset, "mollweide", "originMeanLong", "originMeanLat", verbose = verbose)
  
  dest_subset <- sub_dt[, .(get(taxon), subRegionLong, subRegionLat, subRegionName)]
  setnames(dest_subset, "V1", taxon)
  dest_subset <- unique(dest_subset, by = c(taxon, "subRegionLong", "subRegionLat"))
  dest_points <- get_con_points(dest_subset, "mollweide", "subRegionLong", "subRegionLat", verbose = verbose)
  
  origin_points <- origin_points[!is.na(origin_points)]
  dest_points <- dest_points[!is.na(dest_points)]
  
  catn("Converting points back to data tables.")
  origin <- terra::as.data.frame(origin_points, geom = "xy")
  origin <- as.data.table(origin)
  setnames(origin, "x", "originX")
  setnames(origin, "y", "originY")
  print(length(unique(origin[[taxon]])))
  origin <- unique(origin, by = c(taxon, "originX", "originY"))
  
  dest <- terra::as.data.frame(dest_points, geom = "xy")
  dest <- as.data.table(dest)
  setnames(dest, "x", "destX")
  setnames(dest, "y", "destY")
  dest <- unique(dest, by = c(taxon, "destX", "destY"))
  print(length(unique(dest[[taxon]])))
  
  catn("Merging point data tables.")
  merged_dt <- merge(dest, origin, by = taxon, all = TRUE, allow.cartesian = TRUE)
  merged_dt <- merged_dt[!is.na(destX)]
  merged_dt <- unique(merged_dt, by = c(taxon, "subRegionName"))
  
  vebcat("Unique Taxon:", highcat(length(unique(merged_dt[[taxon]]))))
  
  catn("Creating connections plot.")
  
  min_lim <- 0
  max_lim <- max(merged_dt$connections)
  if (taxon == "species") {
    vis_breaks <- c(min_lim, max_lim)
  } else {
    vis_breaks <- seq(min_lim, max_lim, length.out = 8)
  }
  
  
  fig5 <- ggplot(data = merged_dt) +
    geom_spatvector(data = wm) +
    geom_point(aes(x = originX, y = originY, color = "Origin Country")) +
    geom_point(aes(x = destX, y = destY, color = paste("Arctic", subregion.name))) +
    scale_color_discrete(guide = guide_legend(title.position = "top", ncol = 1)) +
    labs(color = "Points") +
    theme(
      axis.text = element_text(size = 10),
      plot.title = element_text(vjust = 0.5, hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10, hjust = 0.5),
      legend.position = "bottom",
    ) +
    new_scale_color() +
    geom_segment(aes(x = originX, y = originY, xend = destX, yend = destY, color = connections)) +
    ggplot.filler(
      gradient = gradient.config$vis.gradient,
      scale.type = "color-c",
      guide = gradient.config$guide,
      limits = c(min_lim, max_lim),
      breaks = vis_breaks,
      end = ifelse(taxon == "species", 0.5, 1),
      labels = function(x) sprintf("%.0f", round(x, 0)),
      na.value = gradient.config$na.value
    ) +
    labs(
      x = "Longitude",
      y = "Latitude",
      title = if (vis.title) paste0("Potential New Alien ", paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon))), " Richness from origin Country to Arctic Florsitic Province"),
      color = paste0(toupper(substr(taxon, 1, 1)), substr(taxon, 2, nchar(taxon)), " Connections")
    ) +
    theme_minimal() +
    theme.config + theme(
      plot.title = element_text(size = 14),
    )
    
  if (plot.show) print(fig5)
  
  save_ggplot(
    save.plot = fig5, 
    save.name = save.name,
    save.width = 3200, 
    save.height = 2000,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  vebcat("Connection map successfully visualized", color = "funSuccess")
}

visualize_lat_distribution <- function(input.dt, model.scale = "", region.name, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-6", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Absolute Median Latitude plot", color = "funInit")
  
  dt <- copy(input.dt)
  
  lat <- switch(model.scale,
    "log" = log(dt$medianLat),
    "log10" = log10(dt$medianLat),
    "sqrt" = sqrt(dt$medianLat),
    dt$medianLat
  )
  
  overlap <- switch(model.scale,
    "log" = log(dt$overlapRegion),
    "log10" = log10(dt$overlapRegion),
    "sqrt" = sqrt(dt$overlapRegion),
    dt$overlapRegion
  )
  
  min_lim <- 0
  max_lim <- max(dt$medianLat)
  saved_config <- gradient.config$guide
  gradient.config$guide$params$nrow <- NULL
  gradient.config$guide$params$ncol <- 1
  
  fig6 <- ggplot(data = dt, aes(x = lat, y = overlap)) +
    geom_point(aes(color = group)) +
    geom_smooth(method=lm , color="grey", se=TRUE, formula = y ~ x) +
    ggplot.filler(
      gradient = gradient.config$vis.gradient,
      scale.type = "color-d",
      guide = gradient.config$guide,
      begin = 0.5,
      end = 0
    ) +
    labs(
      x = paste0("Absolute Median Latitude", if (model.scale != "") {paste0(" (", model.scale, ")")}),
      y = paste0("Potential Climatic Overlap", if (model.scale != "") {paste0(" (", model.scale, ")")}),
      color = "Taxonomic Group",
      title = if (vis.title) "Influence of Species Latitudinal Ranges on Potential Climatic Overlap"
    ) +
    theme_minimal() +
    theme.config +
    theme(
      legend.position = "right"
    )
    
    if (plot.show) print(fig6)

  if (model.scale != "") save.name <- paste0(save.name, "-", model.scale)
  
  save_ggplot(
    save.plot = fig6, 
    save.name = save.name,
    save.width = 3200, 
    save.height = 2000,
    save.dir = save.dir, 
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title, 
    plot.show = plot.show, 
    verbose = verbose
  )
  
  # reset config
  gradient.config$guide <- saved_config
  
  vebcat("Absolute Median Latitude plot Visualized Successfully", color = "funSuccess")
}


