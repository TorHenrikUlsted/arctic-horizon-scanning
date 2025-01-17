create_custom_labeller <- function(pattern = "_", width = 10) {
  function(labels) {
    wrapped_labels <- lapply(labels, function(col) {
      sapply(col, function(label) {
        wrapped <- paste(strwrap(gsub(pattern, " ", label), width = width), collapse = "\n")
      })
    })
    wrapped_labels
  }
}

#------------------------#
####    Frequency     ####
#------------------------#

visualize_freqpoly <- function(spec.cells, region, region.name, vis.x, vis.color, vis.shade, vis.shade.name, vis.gradient = "viridis-B", vis.title = FALSE, vis.binwidth = 1, vis.x.scale = "log", vis.peak.threshold = NULL, save.dir, save.device = "jpeg", save.unit = "px", plot.show = FALSE, verbose = FALSE) {
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
  # spec.cells <- spec.cells[vis_x != 0, ]

  vis_x <- spec.cells[[vis.x]]

  min_lim <- 1
  max_lim <- max(vis_x, na.rm = TRUE)
  vis_breaks <- c(min_lim, 5, 10, 50, 100, 500, 1000, max_lim)

  if (!is.null(vis.x.scale)) {
    if (vis.x.scale == "log") {
      vis.binwidth <- (vis.binwidth / 10)
    } else if (vis.x.scale == "log10") {
      vis.binwidth <- vis.binwidth / 20
    } else if (vis.x.scale == "sqrt") {
      vis.binwidth <- vis.binwidth / 4
      vis.peak.threshold <- 0.0002
    }
  }
  catn("Using Binwidth:", vis.binwidth)

  # spec.cells, aes(x = vis_x)
  # After_Stat accesses ggplot processed data, count represents the count of data points in each bin of the histogram. Count is cell with that value. So we divide cell by total cells
  fig1A <- ggplot(spec.cells, aes(x = vis_x)) +
    geom_freqpoly(
      binwidth = vis.binwidth,
      aes(y = after_stat(count) / sum(after_stat(count)), color = after_stat(x)),
      na.rm = TRUE,
      show.legend = TRUE
    ) +
    ggplot.filler(
      gradient = config$ggplot$gradient$vis.gradient,
      scale.type = "color-b",
      guide = config$ggplot$gradient$guide,
      breaks = vis_breaks,
      labels = function(x) sprintf("%.0f", round(x, 0)),
      na.value = config$ggplot$gradient$na.value
    ) +
    labs(
      x = paste0("Potential Species Richness ", if (!is.null(vis.x.scale)) {
        paste0("(", vis.x.scale, ")")
      }),
      y = paste0("Proportion of cells in ", region.name),
      title = if (vis.title) paste0("Potential New Alien Species Richness in ", region.name),
      color = "Potential Species Richness"
    ) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.01)) +
    theme_minimal() +
    ggplot.theme() +
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

  library(data.table)

  setDT(peak_data) # Convert to data.table if it's not already

  if (is.null(vis.x.scale)) {
    peak_data[, `:=`(
      cellCount = count,
      count = x
    )]
  } else {
    peak_data[, `:=`(
      cellCount = count,
      count = switch(vis.x.scale,
        "log10" = 10^x,
        "sqrt" = x^2,
        "log" = exp(x),
        x
      )
    )]
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

  #------------------------#
  ####     Freq 1B      ####
  #------------------------#

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

  fig1B <- ggplot(spec.cells, aes(x = vis_x, color = plcol, alpha = plshd)) +
    geom_freqpoly(
      binwidth = 0.1,
      aes(y = after_stat(count) / sum(after_stat(count)))
    ) +
    ggplot.filler(
      gradient = config$ggplot$gradient$vis.gradient,
      guide = config$ggplot$gradient$guide,
      scale.type = "color-d",
      breaks = vis_breaks,
      labels = function(x) sprintf("%.0f", round(x, 0)),
      na.value = config$ggplot$gradient$na.value
    ) +
    scale_alpha_discrete(guide = "none", range = c(vis_amin, vis_amax)) +
    labs(
      x = paste0("Potential Species Richness ", if (!is.null(vis.x.scale)) {
        paste0("(", vis.x.scale, ")")
      }),
      y = paste0("Proportion of cells in the different regions of the ", region.name),
      title = if (vis.title) paste0("Potential New Alien Species Richness in Different Regions of the ", region.name),
      color = "Country",
      alpha = vis.shade.name
    ) +
    scale_y_continuous(limits = c(0, NA), labels = function(x) format(x, big.mark = ",", scientific = FALSE), breaks = seq(0, 1, by = 0.01)) +
    theme_minimal() +
    ggplot.theme() +
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

#------------------------#
####    Hotspots      ####
#------------------------#

visualize_hotspots <- function(raster, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-2", save.device = "jpeg", save.unit = "px", plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Potential Species hotspots", color = "funInit")

  if (save.name != "figure-2") {
    name <- str_to_title(strsplit(save.name, "-")[[1]][[3]])
  } else {
    name <- "Species"
  }

  catn("Getting min and max values.")
  min_lim <- 0
  max_lim <- where.max(raster)
  max_lim <- max(max_lim[, 3])
  vis_breaks <- c(1, 5, 10, 50, 100, 200, 400, 700, 900, max_lim)
  catn("Plotting hotspots.")

  fig2 <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = raster) +
    ggplot.filler(
      gradient = config$ggplot$gradient$vis.gradient,
      guide = config$ggplot$gradient$guide,
      scale.type = "fill-b",
      limits = c(min_lim, max_lim),
      breaks = vis_breaks,
      labels = function(x) sprintf("%.0f", round(x, 0)),
      na.value = config$ggplot$gradient$na.value
    ) +
    labs(
      title = if (vis.title) paste0("Potential New Alien Hotspots for ", name, " in the ", region.name),
      fill = "Potential Species Richness"
    ) +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax),
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() +
    ggplot.theme()

  save_ggplot(
    save.plot = fig2,
    save.name = save.name,
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

  return(fig2)
}

#------------------------#
#### Area of occupany ####
#------------------------#

visualize_paoo <- function(rast, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.row = NULL, vis.col = NULL, vis.title = FALSE, save.dir, save.name = "figure-3", save.device = "jpeg", save.unit = "px", return = FALSE, plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing species with highest Potential Area of occupancy", color = "funInit")

  if (save.name != "figure-3") {
    name <- str_to_title(strsplit(save.name, "-")[[1]][[3]])
  } else {
    name <- "Species"
  }

  catn("Checking crs for raster.")
  rast <- check_crs(rast, projection, projection.method = "near", verbose = verbose)
  catn("Checking crs for region")
  region <- check_crs(region, projection, projection.method = "near", verbose = verbose)

  catn("Plotting hotspots.")

  fig3A <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = rast) +
    ggplot.filler(
      gradient = config$ggplot$gradient$vis.gradient,
      guide = config$ggplot$gradient$guide,
      scale.type = "fill-c",
      breaks = c(0, 1),
      labels = c("0", "1"),
      end = 0.6,
      na.value = config$ggplot$gradient$na.value
    ) +
    labs(
      title = if (vis.title) paste0("Potential New Aliens with Highest Potential\nArea of Occupancy for ", name, " in ", region.name),
      fill = "Potential Area of occupancy"
    ) +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax),
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() +
    ggplot.theme() +
    theme(
      strip.text = element_text(size = 16, face = "italic", lineheight = 0.8),
      panel.spacing = unit(1, "lines")
    )

  custom_labeller <- create_custom_labeller(pattern = config$species$file_separator, width = 10)

  if (!is.null(vis.row)) {
    fig3A <- fig3A +
      facet_wrap(~lyr, nrow = vis.row, ncol = vis.col, labeller = labeller(lyr = custom_labeller))
  }

  if (plot.save) {
    save_ggplot(
      save.plot = fig3A,
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

  vebcat("Successfully visualized highest Potential Area of Occupancy", color = "funSuccess")

  if (return) {
    return(fig3A)
  }
}

#------------------------#
####   Suitability    ####
#------------------------#

visualize_suitability <- function(stack, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.row = NULL, vis.col = NULL, vis.title = FALSE, save.dir, save.name = "figure-3B", save.device = "jpeg", save.unit = "px", plot.save = TRUE, return = FALSE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing suitability plot", color = "funInit")

  if (save.name != "figure-3B") {
    name <- str_to_title(strsplit(save.name, "-")[[1]][[3]])
  } else {
    name <- "Species"
  }

  stack <- check_crs(stack, projection, "bilinear", verbose = verbose)
  region <- check_crs(region, projection, "near", verbose = verbose)

  catn("Acquiring min and max values.")
  # min_lim <- 0.00
  # max_lim <- where.max(stack[[1]])[1, ][[3]] # Get the first row and last item "value"

  max_lim <- max(global(stack, "max", na.rm = TRUE))
  min_lim <- min(global(stack, "min", na.rm = TRUE))

  region_ext <- ext(region)
  vis_breaks <- seq(min_lim, max_lim, length.out = 8)

  catn("Generating plot.")
  fig3B <- ggplot() +
    geom_spatvector(data = region) +
    geom_spatraster(data = stack) +
    ggplot.filler(
      gradient = config$ggplot$gradient$vis.gradient,
      guide = config$ggplot$gradient$guide,
      scale.type = "fill-c",
      limits = c(min_lim, max_lim),
      breaks = vis_breaks,
      labels = function(x) sprintf("%.2f", round(x, 2)),
      na.value = config$ggplot$gradient$na.value
    ) +
    labs(
      title = if (vis.title) paste0("Potential New Alien Climatic suitability\n for ", name, " in ", region.name),
      fill = paste0("Potential Climatic Suitability"),
      color = "Country"
    ) +
    coord_sf(
      xlim = c(extent$xmin, extent$xmax),
      ylim = c(extent$ymin, extent$ymax)
    ) +
    theme_minimal() +
    ggplot.theme() +
    theme(
      strip.text = element_text(size = 16, face = "italic")
    )

  custom_labeller <- create_custom_labeller(pattern = config$species$file_separator, width = 10)

  if (!is.null(vis.row) | !is.null(vis.col)) {
    fig3B <- fig3B + facet_wrap(~lyr, nrow = vis.row, ncol = vis.col, labeller = labeller(lyr = custom_labeller))
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

  if (return) {
    return(fig3B)
  }
}

#----------------------------------#
#### Distribution + suitability ####
#----------------------------------#

### For the figure with paoo and suit - remove space between cols

visualize_dist_suit <- function(stack.paoo, stack.suitability, region, region.name, extent, projection, vis.gradient = "viridis-b", vis.unit = NULL, vis.row, vis.col, vis.title = FALSE, save.dir, save.name = "figure-3D", save.device = "jpeg", save.unit = "px", return = TRUE, plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  paoo_raster <- stack.paoo[[1:vis.row]]
  suit_raster <- stack.suitability[[1:vis.row]]

  fig_paoo <- visualize_paoo(
    rast = paoo_raster,
    region = region,
    region.name = region.name,
    extent = extent,
    projection = projection,
    vis.row = vis.row,
    vis.col = vis.col,
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
    vis.row = vis.row,
    vis.col = vis.col,
    save.dir = save.dir,
    save.device = save.device,
    return = return,
    plot.save = FALSE,
    verbose = verbose
  )

  fig_suit <- fig_suit + ggtitle("Potential Suitability")

  catn("Arranging plots.")

  # margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")

  fig_paoo <- fig_paoo +
    theme(
      plot.margin = margin(0, 0, 0, 0, "cm")
    )

  fig_suit <- fig_suit +
    theme(
      plot.margin = margin(0, 0, 0, 0, "cm"),
    )

  fig3D <- grid.arrange(
    fig_paoo, fig_suit,
    ncol = 2,
    nrow = 1,
    widths = unit(c(0.5, 0.5), "npc")
  )

  h <- 3000 * (5 / 2)
  w <- h * 0.55

  save_ggplot(
    save.plot = fig3D,
    save.name = save.name,
    save.width = w,
    save.height = h,
    save.dir = save.dir,
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title,
    plot.show = plot.show,
    verbose = verbose
  )
}

#------------------------#
####   Composition    ####
#------------------------#

visualize_composition <- function(dt, region.name, vis.x, vis.x.sort, vis.y, vis.fill, vis.group, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.device = "jpeg", save.unit = "px", plot.show = F, verbose = F) {
  vebcat("Visualizing composition plot", color = "funInit")

  dt_copy <- copy(dt)

  vebcat("Number of Sub Regions:", highcat(length(unique(dt_copy[[vis.x]]))))

  # Order dt by vis.x.sort column
  setorderv(dt_copy, vis.x.sort)

  # Get ordered levels for x-axis
  ordered_levels <- unique(dt_copy[[vis.x]])

  # Convert x variable to ordered factor
  dt_copy[, (vis.x) := factor(get(vis.x), levels = ordered_levels)]

  # Convert fill variable to factor and order if needed
  dt_copy[, (vis.fill) := factor(get(vis.fill))]

  vebprint(levels(dt_copy[[vis.fill]]), verbose, "Initial fill levels:")

  # Debug output
  vebprint(
    head(dt_copy[, .(
      subRegion = get(vis.x),
      richness = get(vis.y),
      taxon = get(vis.fill)
    )]),
    verbose, "Sample of data:"
  )

  # Apply any additional grouping
  dt_copy <- get_order_group(dt_copy, verbose = verbose)

  vebprint(levels(dt_copy[[vis.fill]]), verbose, "Final fill levels:")

  # Create base plot ensuring proper variable mapping
  p1 <- geom_bar(stat = "identity", position = "stack")

  # Save config
  saved_config_params <- config$ggplot$gradient$guide$params
  config$ggplot$gradient$guide$params$nrow <- NULL
  config$ggplot$gradient$guide$params$ncol <- 1

  p2 <- ggplot.filler(
    gradient = config$ggplot$gradient$vis.gradient,
    guide = list(ncol = 2),
    scale.type = "fill-d",
    na.value = config$ggplot$gradient$na.value
  )
  p4 <- labs(
    x = "Floristic Province",
    y = paste0(
      "Potential Relative ",
      paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill))),
      " Richness"
    ),
    title = if (vis.title) {
      paste0(
        "Potential New Alien ",
        paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill))),
        " Composition in ",
        region.name
      )
    },
    fill = paste0(toupper(substr(vis.fill, 1, 1)), substr(vis.fill, 2, nchar(vis.fill)))
  )
  p5 <- theme_minimal()
  p6 <- theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 16),
    legend.position = "right"
  )

  fig4A <- ggplot(dt_copy, aes(
    x = .data[[vis.x]],
    y = .data[[vis.y]],
    fill = .data[[vis.fill]]
  )) +
    p1 +
    p2 +
    p4 +
    p5 +
    ggplot.theme() +
    p6

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
      gradient = config$ggplot$gradient$vis.gradient,
      scale.type = "fill-d",
      guide = guide_legend(ncol = 1),
      na.value = config$ggplot$gradient$na.value
    ) +
    p4 +
    p5 +
    ggplot.theme() +
    p6

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
  config$ggplot$gradient$guide$params <- saved_config_params

  vebcat("Composition plot successfully visualized", color = "funSuccess")
}

#------------------------#
####   Connections    ####
#------------------------#

visualize_connections <- function(dt, taxon, centroid = FALSE, region.name, subregion.name, projection, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-5", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Connections map", color = "funInit")

  if (taxon == "species") taxon <- "cleanName"

  wm <- get_world_map(projection = projection, verbose = verbose)
  # wm <- wm[wm$level3Name != "Antarctica"]

  sub_dt <- copy(dt)

  if (taxon == "cleanName") {
    origin_subset <- sub_dt[, .(cleanName, originLong, originLat, connections, originCountry, level2Code)]
  } else {
    origin_subset <- sub_dt[, .(cleanName, get(taxon), originLong, originLat, connections, originCountry, level2Code)]
    setnames(origin_subset, "V2", taxon)
  }
  origin_subset <- unique(origin_subset, by = c("cleanName", "originCountry"))

  if (centroid) {
    origin_subset[, level3Name := originCountry]
    print(origin_subset)
    wgsrpd_centroids <- get_wgsrpd_polygon_centroid(origin_subset)
    print(wgsrpd_centroids)
    origin_subset <- origin_subset[wgsrpd_centroids, on = "level3Name"]
    print(origin_subset)
    data.table::setnames(origin_subset, c("originLong", "originLat"), c("level3Long", "level3Lat"))
    data.table::setnames(origin_subset, c("centroidLong", "centroidLat"), c("originLong", "originLat"))
  }

  origin_points <- get_con_points(origin_subset, projection, "originLong", "originLat", verbose = verbose)

  if (taxon == "cleanName") {
    dest_subset <- sub_dt[, .(cleanName, subRegionLong, subRegionLat, subRegionName)]
  } else {
    dest_subset <- sub_dt[, .(cleanName, get(taxon), subRegionLong, subRegionLat, subRegionName)]
    setnames(dest_subset, "V2", taxon)
  }
  dest_subset <- unique(dest_subset, by = c("cleanName", "subRegionName"))
  dest_points <- get_con_points(dest_subset, projection, "subRegionLong", "subRegionLat", verbose = verbose)

  # Remove NA
  origin_points <- origin_points[!is.na(origin_points)]
  dest_points <- dest_points[!is.na(dest_points)]

  catn("Converting points back to data tables.")
  origin <- terra::as.data.frame(origin_points, geom = "xy")
  origin <- as.data.table(origin)
  setnames(origin, "x", "originX")
  setnames(origin, "y", "originY")
  origin <- unique(origin, by = c("cleanName", "originX", "originY"))

  dest <- terra::as.data.frame(dest_points, geom = "xy")
  dest <- as.data.table(dest)
  setnames(dest, "x", "destX")
  setnames(dest, "y", "destY")
  dest <- unique(dest, by = c("cleanName", "destX", "destY"))

  # origin <- unique(origin, by = c("cleanName", "originX", "originY"))
  # dest <- unique(dest, by = c("cleanName", "destX", "destY"))
  origin <- unique(origin, by = c("originX", "originY"))
  dest <- unique(dest, by = c("destX", "destY"))

  # catn("Merging point data tables.")
  # merged_dt <- dest[origin, on = "cleanName", allow.cartesian = TRUE, nomatch = NULL]
  # merged_dt <- merged_dt[!is.na(destX)]
  #
  # catn("Calculating number of connections")
  # merged_dt[, connections := uniqueN(cleanName), by = c(taxon, "subRegionName", "originCountry")]
  # merged_dt <- unique(merged_dt, by = c(taxon, "subRegionName", "originCountry"))
  #
  # vebcat("Unique Taxon:", highcat(length(unique(merged_dt[[taxon]]))))
  #
  # catn("Creating connections plot.")

  # min_lim <- 0
  # max_lim <- max(merged_dt$connections)
  # if (taxon == "cleanName") {
  #   vis_breaks <- c(min_lim, max_lim)
  # } else {
  #   vis_breaks <- seq(min_lim, max_lim, length.out = 8)
  # }

  if (taxon == "cleanName") taxon <- "species"
  region_name <- paste("Arctic", subregion.name)
  fig5 <- ggplot() + # data = merged_dt
    geom_spatvector(data = wm) +
    # geom_segment(aes(x = originX, y = originY,
    #                  xend = destX, yend = destY,
    #                  color = connections),
    #              alpha = if(taxon == "species") 0.1 else 1, # Set a base transparency for all lines
    #              linewidth = 0.4) +
    # ggplot.filler(
    #   gradient = config$ggplot$gradient$vis.gradient,
    #   scale.type = "color-c",
    #   guide = config$ggplot$gradient$guide,
    #   limits = if(taxon == "species") c(1,1) else c(min_lim, max_lim),
    #   breaks = if(taxon == "species") c(1) else vis_breaks,
    #   end = if(taxon == "species") 0 else 1,
    #   labels = if(taxon == "species") c("1") else function(x) sprintf("%.0f", round(x, 0)),
    #   na.value = config$ggplot$gradient$na.value
    # ) +
    #
    # theme(
    #   plot.title = element_text(size = 14),
    # ) +
    # new_scale_color() +
    geom_point(data = origin, aes(x = originX, y = originY, color = "Species Occurrence Mean Range"), size = 0.05) +
    geom_point(data = dest, aes(x = destX, y = destY, color = "Arctic Floristic Province"), size = 0.5) +
    geom_spatvector(data = wm, fill = NA, color = "black", linewidth = 0.1) + # Add country outlines
    scale_color_manual(
      values = c(
        "Species Occurrence Mean Range" = "darkgreen",
        "Arctic Floristic Province" = "#E66101"
      ),
      guide = guide_legend(title.position = "top", ncol = 1)
    ) +
    labs(
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    theme_minimal() +
    ggplot.theme() +
    theme(
      axis.text = element_text(size = 0),
      plot.title = element_text(vjust = 0, hjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 0, hjust = 0.5),
      legend.key.height = unit(0, "cm"),
      legend.position = "bottom",
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.spacing.y = unit(0, "pt"),
      legend.spacing.x = unit(0, "pt"),
      legend.box.spacing = unit(0, "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt")
    )

  if (plot.show) print(fig5)

  save_ggplot(
    save.plot = fig5,
    save.name = save.name,
    save.width = 2000,
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

visualize_spec_ranges <- function(dt, taxon, centroid = FALSE, multiple = FALSE, region, region.name, subregion.name, projection, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-5", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Connections map", color = "funInit")

  if (taxon == "species") taxon <- "cleanName"

  wm <- get_world_map(projection = projection, verbose = verbose)
  region <- load_region(region)
  region <- region[, c("floregName", "floregLong", "floregLat")]
  region <- check_crs(region, projection)
  # wm <- wm[wm$level3Name != "Antarctica"]
                          
  if (multiple) subregion <- region[region$floregName == subregion.name]
  
  catn("Number of rows in input data to function:", highcat(nrow(dt)))
  if (multiple) catn("Any not the input subregion:", highcat(any(dt$subRegionName != subregion.name)))
  
  converted_points <- convert_con_points(
    dt = dt,
    taxon = taxon,
    projection = projection,
    centroid = centroid,
    verbose = verbose
  )
  
  if (multiple) catn("Any not the input subregion after conversion:", highcat(any(dt$subRegionName != subregion.name)))
  
  catn("Number of rows in converted data:", highcat(nrow(converted_points$origin)))

  if (taxon == "cleanName") taxon <- "species"

  # Create base map
  fig5 <- ggplot() +
    geom_spatvector(data = wm) +
    geom_point(data = converted_points$origin, aes(x = originX, y = originY, color = "Species Points"), size = 0.35, stroke = 0, shape = 16) +
    geom_spatvector(
      data = region,
      aes(fill = "Subregions"),
      color = "transparent"
    ) +
    geom_spatvector(data = wm, fill = NA, color = "black", linewidth = 0.1)
  
  # Create the zoomed plot of the subregion
  if (multiple) {
    
    zoom_grob <- create_zoom_annotation(
      basemap = wm,
      region = subregion,
      points = converted_points$origin
    )
    
    # Modified positioning for the inset
    wm_ext <- terra::ext(wm)
    region_ext <- terra::ext(subregion)
    
    inset_pos <- position_zoom_annotation(
      top = -0.05, 
      right = -0.1, 
      plot.ext = wm_ext,
      region.ext = zoom_grob$extents,
      width = 0.35, 
      height = 0.25
    )
    # xmin & ymin = bottom-left
    segment_extents <- data.table(
      x = c(region_ext[1], region_ext[2]), # left-top, right-bottom
      y = c(region_ext[4], region_ext[3]),
      xend = c(inset_pos$xmin, inset_pos$xmax), # left-bottom, right-bottom
      yend = c(inset_pos$ymin, inset_pos$ymin)
    )
    
    fig5 <- fig5 +
      geom_spatvector(
        data = subregion,
        aes(fill = "Selected Subregion"),
        color = "transparent"
      ) +
      geom_segment(
        data = segment_extents,
        aes(x = x, y = y, xend = xend, yend = yend),
        linetype = "dashed",
        color = "gray10",
        linewidth = 0.3
      ) +
      annotation_custom(
        grob = zoom_grob$grob,
        xmin = inset_pos$xmin,
        xmax = inset_pos$xmax,
        ymin = inset_pos$ymin,
        ymax = inset_pos$ymax
      )
  }

  # Add points and outlines on top of the basemaps
  fig5 <- fig5 +
    scale_color_manual(
      values = c(
        "Species Points" = "darkgreen"
      ),
      labels = c(
        "Species Points" = "Species Centroids within Botanical Countries"
      ),
      guide = guide_legend(
        title.position = "top", ncol = 1, 
        order = 1, override.aes = list(size = 3)
      )
    )
  
    # Add the appropriate fill scale based on whether multiple is true
    if (multiple) {
      fig5 <- fig5 +
        scale_fill_manual(
          values = c(
            "Subregions" = "#4682B4",
            "Selected Subregion" = "#0013ff"
          ),
          labels = c(
            "Subregions" = "Circumpolar Arctic",
            "Selected Subregion" = subregion.name
          ),
          guide = guide_legend(title.position = "top", ncol = 1, order = 2)
        )
    } else {
      fig5 <- fig5 +
        scale_fill_manual(
          values = c(
            "Subregions" = "#4682B4"
          ),
          labels = c(
            "Circumpolar Arctic"
          ),
          guide = guide_legend(title.position = "top", ncol = 1, order = 2)
        )
    }

  # Add the theme to the plot
  fig5 <- fig5 +
    labs(
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    theme_minimal() +
    ggplot.theme() +
    theme(
      axis.text = element_text(size = 0),
      plot.title = element_text(vjust = 0, hjust = 0.5),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.title = element_text(size = 0, hjust = 0.5),
      legend.text = element_text(size = 8),
      legend.key.height = unit(0.5, "cm"),
      legend.key.width = unit(0.5, "cm"),
      legend.spacing.y = unit(0, "pt"),
      legend.spacing.x = unit(0, "pt"),
      legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.spacing = unit(0, "pt")
    )

  if (plot.show) print(fig5)

  save_ggplot(
    save.plot = fig5,
    save.name = save.name,
    save.width = 2000,
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

visualize_spec_ranges_by_subregion <- function(dt, taxon, centroid = FALSE, multiple = FALSE, region, region.name, subregion.name, projection, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-5", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  save.dir <- file.path(save.dir, "connections")

  # Get unique subregions
  subregions <- unique(dt$subRegionName)

  # First create the overall plot with all subregions
  visualize_spec_ranges(
    dt = dt,
    taxon = taxon,
    centroid = centroid,
    multiple = FALSE,
    region = region,
    region.name = region.name,
    subregion.name = subregion.name,
    projection = projection,
    vis.gradient = vis.gradient,
    vis.title = vis.title,
    save.dir = save.dir,
    save.name = paste0(save.name, "-all"),
    save.device = save.device,
    save.unit = save.unit,
    plot.save = plot.save,
    plot.show = plot.show,
    verbose = verbose
  )
  
  catn("Number of total rows of input data:", highcat(nrow(dt)))
  vebprint(names(dt), text = "Names of input data:")
  catn("Number of subregions:", highcat(length(subregions)))
  
  if (multiple) {
    # Then create individual plots for each subregion
    for (i in 1:length(subregions)) {
      subregion <- subregions[[i]]
      # Filter data for current subregion
      subregion_dt <- dt[subRegionName == subregion]
      
      catn("Subsetting data for:", highcat(subregion))
      catn("Number of rows in subregion:", highcat(nrow(subregion_dt)))
      catn("Proportion of rows subregion/total rows:", highcat(nrow(subregion_dt) / nrow(dt)))
      
      # Skip if no data for this subregion
      if (nrow(subregion_dt) == 0) {
        warning(paste("No data for subregion:", subregion))
        next
      }
      
      # Create plot for this subregion
      visualize_spec_ranges(
        dt = subregion_dt,
        taxon = taxon,
        centroid = centroid,
        multiple = TRUE,
        region = region,
        region.name = region.name,
        subregion.name = subregion,
        projection = projection,
        vis.gradient = vis.gradient,
        vis.title = vis.title,
        save.dir = save.dir,
        save.name = paste0(save.name, "-", gsub(" ", "-", subregion)),
        save.device = save.device,
        save.unit = save.unit,
        plot.save = plot.save,
        plot.show = plot.show,
        verbose = verbose
      )
    }
  }
}

#------------------------------------------#
####   Distribution - beta Regression   ####
#------------------------------------------#

visualize_gamlss <- function(dt, model, region.name, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-6", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Absolute Median Latitude plot", color = "funInit")
  data <- copy(dt)

  # Create prediction data table
  lat_seq <- seq(min(data$medianLat), max(data$medianLat), length.out = 100)
  pred_dt <- data.table(medianLat = lat_seq)

  # Get coefficients from model
  nu_coef <- coef(model$models[[1]], what = "nu") # Zero-inflation coefficients
  mu_coef <- coef(model$models[[1]], what = "mu") # Mean coefficients

  # Calculate predicted components
  logit <- function(x) exp(x) / (1 + exp(x))

  # Add all predictions
  pred_dt[, `:=`(
    prob_nonzero = 1 - logit(nu_coef[1] + (nu_coef[2] * medianLat)),
    expected_mean = logit(mu_coef[1] + (mu_coef[2] * medianLat))
  )]
  pred_dt[, combined := prob_nonzero * expected_mean]

  # Melt to long format
  pred_long <- melt(pred_dt,
    id.vars = "medianLat",
    measure.vars = c("prob_nonzero", "expected_mean", "combined"),
    variable.name = "component",
    value.name = "value"
  )
  # Create the plot
  fig6A <- ggplot() +
    # Add raw data points
    geom_point(
      data = data,
      aes(x = medianLat, y = overlapRegion),
      alpha = 0.2, color = "grey50", size = 1
    ) +
    # Add model predictions
    geom_line(
      data = pred_long,
      aes(x = medianLat, y = value, color = component),
      linewidth = 1
    ) +
    # Customize colors and labels
    scale_color_manual(
      values = c(
        "prob_nonzero" = "#8884d8",
        "expected_mean" = "#82ca9d",
        "combined" = "#ff7300"
      ),
      labels = c(
        "prob_nonzero" = "Prob(Overlap > 0)",
        "expected_mean" = "E(Overlap | Overlap > 0)",
        "combined" = "Combined Expectation"
      )
    ) +
    # Customize theme and labels
    theme_bw() +
    labs(
      x = "Species Absolute Median Latitude (°)",
      y = paste0("Proportional Niche Overlap of ", region.name),
      if (vis.title) title <- "Niche Overlap vs Absolute Median Latitude",
      if (vis.title) subtitle <- "GAMLSS model components with occurrence data",
      color = "Component"
    ) +
    theme(legend.position = "bottom")

  if (plot.show) print(fig6A)

  save_ggplot(
    save.plot = fig6A,
    save.name = paste0(save.name, "A"),
    save.width = 3200,
    save.height = 2000,
    save.dir = save.dir,
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title,
    plot.show = plot.show,
    verbose = verbose
  )

  # For the two-panel version:
  fig6B <- ggplot() +
    geom_point(
      data = data,
      aes(x = medianLat, y = overlapRegion),
      alpha = 0.2, color = "grey50", size = 1
    ) +
    geom_line(
      data = pred_dt,
      aes(x = medianLat, y = combined),
      color = "#ff7300", linewidth = 1
    ) +
    theme_bw() +
    labs(
      x = "Species Absolute Median Latitude",
      y = "Proportional Arctic Niche Overlap",
      title = "A. Predicted Overlap with "
    )

  fig6B2 <- ggplot() +
    geom_line(
      data = pred_dt,
      aes(x = medianLat, y = prob_nonzero),
      color = "#8884d8", linewidth = 1
    ) +
    theme_bw() +
    labs(
      x = "Species Absolute Median Latitude",
      y = "Probability of Non-Zero Overlap",
      title = "B. Probability of Any Overlap"
    ) +
    scale_y_continuous(limits = c(0, 1))

  fig6B <- fig6B / fig6B2 + plot_layout(heights = c(3, 2))

  if (plot.show) print(fig6B)

  save_ggplot(
    save.plot = fig6B,
    save.name = paste0(save.name, "B"),
    save.width = 3200,
    save.height = 2000,
    save.dir = save.dir,
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title,
    plot.show = plot.show,
    verbose = verbose
  )

  vebcat("Absolute Median Latitude plot Visualized Successfully", color = "funSuccess")
}

#------------------------------#
####   Distribution - GAM   ####
#------------------------------#
visualize_gam <- function(dt, model, region.name, vis.gradient = "viridis-b", vis.title = FALSE, save.dir, save.name = "figure-7", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing Latitude GAM plot", color = "funInit")
  data <- copy(dt)
  
  # Create prediction data frame
  lat_seq <- seq(min(dt$centroidLatitude), max(dt$centroidLatitude), length.out = 100)
  pred_df <- data.table(centroidLatitude = lat_seq)
  
  # Get predictions with standard errors
  predictions <- predict(model$full, newdata = pred_df, se.fit = TRUE, type = "response")
  
  # Create prediction data table with confidence intervals
  pred_dt <- data.table(
    latitude = lat_seq,
    fit = predictions$fit,
    se = predictions$se.fit
  )
  
  # Add confidence intervals
  pred_dt[, `:=`(
    lower = fit - (1.96 * se),
    upper = fit + (1.96 * se)
  )]
  
  # Create the plot
  fig7 <- ggplot() +
    # Add confidence interval ribbon
    geom_ribbon(
      data = pred_dt,
      aes(x = latitude, ymin = lower, ymax = upper),
      alpha = 0.2,
      fill = "#ff7300"
    ) +
    # Add raw data points
    geom_point(
      data = data,
      aes(x = centroidLatitude, y = overlapRegion),
      alpha = 0.2,
      color = "grey50",
      size = 1
    ) +
    # Add fitted line
    geom_line(
      data = pred_dt,
      aes(x = latitude, y = fit),
      color = "#ff7300",
      linewidth = 1
    ) +
    # Add reference line at equator
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey50",
      alpha = 0.5
    ) +
    theme_bw() +
    labs(
      x = "Species Centroid Latitude (°)",
      y = paste0("Species Proportional Niche Overlap with ", region.name),
      title = if(vis.title) "Niche Overlap vs Latitude" else NULL,
      subtitle = if(vis.title) "GAM with beta regression" else NULL
    )
  
  if (plot.show) print(fig7)
  
  save_ggplot(
    save.plot = fig7,
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
  
  vebcat("Latitude GAM plot Visualized Successfully", color = "funSuccess")
}


#------------------------#
####      Sankey      ####
#------------------------#

visualize_sankey <- function(dt, taxon, vis.gradient = "viridis-b", vis.title = FALSE, region.name = "Region", subregion.name = "Sub Region", save.dir, save.name = "figure-8", save.device = "jpeg", save.unit = "px", plot.save = TRUE, plot.show = FALSE, verbose = FALSE) {
  vebcat("Visualizing data in a sankey plot", color = "funInit")

  dt_sank <- copy(dt)

  catn("Number of origin countries:", highcat(length(unique(dt_sank$origin))))
  catn("Number of", paste0(subregion.name, "s:"), highcat(length(unique(dt_sank$destination))))

  # Function to wrap long text
  wrap_text <- function(text, width = 20) {
    sapply(strwrap(text, width = width, simplify = FALSE),
      paste,
      collapse = "\n"
    )
  }

  # Apply text wrapping to long names
  dt_sank$origin_wrapped <- wrap_text(dt_sank$origin)
  dt_sank$destination_wrapped <- wrap_text(dt_sank$destination)

  saved_config <- config$ggplot$gradient$guide
  config$ggplot$gradient$guide$label.position <- "right"
  config$ggplot$gradient$guide$nrow <- NULL
  config$ggplot$gradient$guide$ncol <- 1

  catn("Creating sankey plot.")

  fig8 <- ggplot(
    data = dt_sank,
    aes(
      axis1 = origin_wrapped,
      axis2 = destination_wrapped,
      y = relativeRichness
    )
  ) +
    geom_flow(aes(fill = destination),
      width = 0.3,
      alpha = 0.7,
      curve_type = "quintic"
    ) + # Smoother curves
    geom_stratum(
      width = 0.3,
      fill = "grey95",
      color = "grey40"
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 3,
      hjust = 0.5,
      check_overlap = TRUE
    ) +
    scale_x_discrete(
      limits = c("Origin Country", subregion.name),
      expand = c(0.05, 0.05)
    ) +
    ggplot.filler(
      gradient = config$ggplot$gradient$vis.gradient,
      scale.type = "fill-d",
      guide = config$ggplot$gradient$guide,
    ) +
    labs(
      y = paste0(
        "Relative ", toupper(substr(taxon, 1, 1)),
        substr(taxon, 2, nchar(taxon)), " Richness"
      ),
      if (vis.title) {
        title <- paste0(
          "Potential New Alien Relative ",
          toupper(substr(taxon, 1, 1)),
          substr(taxon, 2, nchar(taxon)),
          " Richness"
        )
      },
      if (vis.title) subtitle <- paste0("From Origin Country to ", subregion.name),
      fill = paste0(subregion.name, "s in ", region.name)
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 8), # Smaller axis text
      axis.text.x = element_text(angle = 0), # Horizontal axis labels
      plot.title = element_text(
        size = 14, face = "bold",
        vjust = 2, hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 12, hjust = 0.5,
        vjust = 2, color = "grey40"
      ),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      panel.grid = element_blank(),
      # Add more space at the bottom for labels
      plot.margin = margin(t = 30, r = 30, b = 50, l = 30)
    )

  config$ggplot$gradient$guide <- saved_config

  if (plot.show) print(fig8)

  save_ggplot(
    save.plot = fig8,
    save.name = save.name,
    save.width = 3840,
    save.height = 3500,
    save.dir = save.dir,
    save.device = save.device,
    save.unit = save.unit,
    vis.title = vis.title,
    plot.show = plot.show,
    verbose = verbose
  )

  vebcat("Sankey plot successfully visualized", color = "funSuccess")
}
