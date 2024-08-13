source_all("./src/visualize/components")

visualize_sequence <- function(out.dir = "./outputs/visualize", res.unknown, res.known, shape, hv.dir, hv.method = "box", vis.projection = "longlat", vis.title = TRUE, vis.region.name = "Region", vis.subregion.name = "Sub Region", vis.composition.taxon = "order", vis.save.device = "jpeg", vis.save.unit = "px", plot.show = FALSE, verbose = FALSE) {
  
  # Test values
  out.dir = "./outputs/visualize/glonaf"
  res.unknown = "glonaf" 
  res.known = "arctic"
  shape = "./outputs/setup/region/cavm-noice/cavm-noice.shp" 
  hv.dir = "./outputs/hypervolume/glonaf" 
  hv.method = "box"
  vis.projection = "laea"
  vis.title = TRUE
  vis.region.name = "Region"
  vis.subregion.name = "Sub Region"
  vis.composition.taxon = "order"
  vis.save.device = "jpeg"
  vis.save.unit = "px"
  plot.show = FALSE
  verbose = FALSE
  
  
  ##########################
  #       Load dirs        #
  ##########################
  hv_proj_dir <- paste0(hv.dir, "/", hv.method, "-sequence/projections")
  vis_dir <- paste0(out.dir, "/", hv.method, "-sequence")
  stats_dir <- paste0(vis_dir, "/stats")
  log_dir <- paste0(vis_dir, "/logs")
  result_dir <- paste0(vis_dir, "/results")
  rast_dir <- paste0(log_dir, "/rasters")
  plot_dir <- paste0(vis_dir, "/plots")

  # create dirs
  create_dir_if(stats_dir, log_dir, result_dir, rast_dir, plot_dir)

  ##########################
  #       Load files       #
  ##########################
  # Set up error handling
  warn_file <- paste0(log_dir, "/", hv.method, "-warning.txt")
  err_file <- paste0(log_dir, "/", hv.method, "-error.txt")
  stats_file <- paste0(hv.dir, "/", hv.method, "-sequence/stats/stats.csv")
  expected_res_file <- paste0("./outputs/filter/", res.unknown, "/", res.unknown, "-absent-final.csv")
  unknown_chunk_dir <- paste0("./outputs/filter/", res.unknown, "/chunk/species")
  known_res_file <- paste0("./outputs/filter/", res.known, "/", res.known, "-present-final.csv")
  res_file <- paste0(hv.dir, "/", hv.method, "-sequence/stats/stats.csv")
  if (vis.title) {
    plot_dir <- paste0(plot_dir, "/title/", vis.save.device)
    existing_plots <- basename(list.files(plot_dir))
  } else {
    plot_dir <- paste0(plot_dir, "/no-title/", vis.save.device)
    existing_plots <- basename(list.files(plot_dir))
  }

  create_file_if(warn_file, err_file)

  ##########################
  #     Load varaibles     #
  ##########################

  sp_dirs <- list.dirs(hv_proj_dir, full.names = TRUE)

  # Remove the directory name
  sp_dirs <- sp_dirs[-1]

  region <- load_region(shape)
  region <- check_crs(region, vis.projection)

  ##########################
  #       Check data       #
  ##########################

  # Clean species that are not supposed to be there
  stats_clean <- clean_output_species(
    out.dir = stats_dir,
    dt.result = res_file,
    dt.expected = expected_res_file,
    clean.dirs = hv_proj_dir,
    verbose = FALSE
  )

  # Get stats csv file
  sp_stats <- filter_stats_data(
    stats_clean,
    known_res_file,
    unknown_chunk_dir,
    stats_dir,
    hv.dir = hv.dir,
    hv.method = hv.method,
    verbose = verbose
  )

  catn("Found", highcat(length(sp_dirs)), "projections and", highcat(length(unique(sp_stats$cleanName))), "Species.")

  if (length(sp_dirs) > length(unique(sp_stats$cleanName))) {
    stop("Found more projections than species.")
  } else if (length(sp_dirs) < length(unique(sp_stats$cleanName))) {
    stop("Found more species than projections.")
  }

  ##########################
  #     Get cell data      #
  ##########################

  template_file <- list.dirs(hv_proj_dir)[2]
  template_file <- paste0(template_file, "/0.5-inclusion.tif")

  # Get which cells the sub regions are in
  region_cell <- get_region_cells(
    shape = shape,
    template.filename = template_file,
    out.dir = result_dir,
    verbose = FALSE
  )
  
  # Split filenames into groups
  sp_group_dirs <- split_spec_by_group(sp_dirs, sp_stats, "cleanName", verbose)
  
  group_richness <- list()
  
  for (group in names(sp_group_dirs)) {
    group_dirs <- sp_group_dirs[[group]]$filename
    # Get the cell richness per group
    inc_dt <- parallel_spec_handler(
      spec.dirs = group_dirs,
      dir = paste0(log_dir, "/", group, "-cell"),
      shape = shape,
      hv.project.method = "0.5-inclusion",
      fun = get_inclusion_cell,
      batch = TRUE,
      node.log.append = FALSE,
      verbose = verbose
    )
    
    # combine cell richness and region cell occupancy
    region_richness_cell <- merge(inc_dt, region_cell, by = "cell", all = TRUE)
    
    group_richness[[group]] <- region_richness_cell
  }
  
  summed_richness <- copy(group_richness[[1]])
  
  for (i in 2:length(group_richness)) {
    group_data <- group_richness[[i]]
    
    summed_richness[group_data, on = "cell", cellRichness := 
                fifelse(is.na(x.cellRichness) & is.na(i.cellRichness), NA_real_,
                        fifelse(is.na(x.cellRichness), i.cellRichness,
                                fifelse(is.na(i.cellRichness), x.cellRichness,
                                        x.cellRichness + i.cellRichness)))
    ]
    
  }
  


  ##########################
  #        Figure 1        #
  ##########################

  fig_name <- paste0("figure-1A", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Frequency Figures.", color = "indicator")
  } else {
    freq_stack <- region_richness_cell[!is.na(region_richness_cell$cellRichness)]
    freq_stack <- freq_stack[!is.na(freq_stack$subRegionName)]

    visualize_freqpoly(
      spec.cells = freq_stack,
      region = region,
      region.name = vis.region.name,
      vis.x = "cellRichness",
      vis.title = vis.title,
      vis.color = "country",
      vis.shade = "subRegionName",
      vis.shade.name = vis.subregion.name,
      vis.binwidth = 1,
      vis.x.scale = "log",
      vis.peak.threshold = 0.001,
      save.dir = plot_dir,
      save.device = vis.save.device,
      verbose = verbose
    )

    mdwrite(
      config$files$post_seq_md,
      text = "2;Potential Frequency",
    )

    mdwrite(
      config$files$post_seq_md,
      text = "See figure 1A descriptive for result numbers."
    )

    rm(freq_stack)
    invisible(gc())
  }

  ##########################
  #        Figure 2        #
  ##########################

  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots

  fig_name <- paste0("figure-2", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Hotspots Figure.", color = "indicator")
  } else {
    # Convert to raster by using a template
    hotspot_raster <- convert_template_raster(
      input.values = region_richness_cell$cellRichness,
      template.filename = template_file,
      projection = out_projection,
      projection.method = "near",
      out.dir = rast_dir,
      verbose = verbose
    )

    world_map <- get_world_map(projection = out_projection)

    visualize_hotspots(
      rast = hotspot_raster,
      region = world_map,
      extent = region_ext,
      region.name = vis.region.name,
      projection = out_projection,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      vis.title = vis.title,
      verbose = verbose
    )

    res_hotspots_nums <- region_richness_cell[, .SD[which.max(cellRichness)], by = subRegionName]
    res_hotspots_nums <- res_hotspots_nums[!is.na(res_hotspots_nums$subRegionName)]
    res_hotspots_nums <- res_hotspots_nums[, .(cellRichness, subRegionName, country)]
    setnames(res_hotspots_nums, "cellRichness", "richnessPeaks")
    setorder(res_hotspots_nums, -richnessPeaks)

    mdwrite(
      config$files$post_seq_md,
      text = "2;Hotspots",
      data = res_hotspots_nums,
    )

    rm(hotspot_raster, world_map, res_hotspots_nums)
    invisible(gc())
  }

  ##########################
  #        Figure 3A       #
  ##########################

  fig_name <- paste0("figure-3A", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Potential Area of Occupancy Figure.", color = "indicator")
  } else {
    paoo_file <- paste0(log_dir, "/area-of-occupancy/0.5-inclusion/area-of-occupancy.csv")

    paoo_files <- parallel_spec_handler(
      spec.dirs = sp_dirs,
      shape = shape,
      dir = dirname(dirname(paoo_file)),
      hv.project.method = "0.5-inclusion",
      col.n = "species-9",
      out.order = "-TPAoO",
      fun = get_paoo,
      verbose = verbose
    )

    paoo <- stack_projections(
      filenames = paoo_files$filename,
      projection = out_projection,
      projection.method = "near",
      out.dir = paste0(log_dir, "/stack-aoo"),
      binary = TRUE,
      verbose = verbose
    )

    world_map <- get_world_map(projection = out_projection)

    visualize_paoo(
      rast = paoo,
      region = world_map,
      region.name = vis.region.name,
      extent = region_ext,
      projection = out_projection,
      vis.wrap = 3,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      vis.title = vis.title,
      verbose = verbose
    )

    paoo_md <- paoo_files[, 1:3]

    mdwrite(
      config$files$post_seq_md,
      text = "2;Highest Total Potential Area of Occupancy",
      data = paoo_md,
    )


    rm(paoo_files, paoo, paoo_md, world_map)
    invisible(gc())
  }

  ##########################
  #        Figure 3B       #
  ##########################

  fig_name <- paste0("figure-3B", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Suitability Figure.", color = "indicator")
  } else {
    # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
    # Get the max values of all probability raster files
    prob_mean <- parallel_spec_handler(
      spec.dirs = sp_dirs,
      shape = shape,
      dir = paste0(log_dir, "/suitability"),
      hv.project.method = "probability",
      col.n = "species-9",
      out.order = "-totalMean",
      fun = get_prob_stats,
      verbose = verbose
    )

    # Read the rasters with the highest max values
    prob_stack <- stack_projections(
      filenames = prob_mean$filename,
      projection = out_projection,
      projection.method = "bilinear",
      out.dir = paste0(log_dir, "/stack-suitability-", "totalMean"),
      verbose = verbose
    )

    world_map <- get_world_map(projection = out_projection)

    visualize_suitability(
      stack = prob_stack,
      region = world_map,
      region.name = vis.region.name,
      extent = region_ext,
      projection = out_projection,
      vis.unit = "mean",
      vis.wrap = 3,
      vis.title = vis.title,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      verbose = verbose
    )

    prob_vals_md <- copy(prob_mean)
    prob_vals_md <- prob_vals_md[, filename := NULL]

    prob_vals_md <- prob_vals_md[, 1:5]

    mdwrite(
      config$files$post_seq_md,
      text = paste0(
        "2;Highest Potential Climatic Suitability\n\n",
        "See figure 3D for comparison between **potential area of occupancy** and **potential climatic suitability**  "
      ),
      data = prob_vals_md
    )

    rm(prob_mean, prob_stack, world_map, prob_vals_md)
    invisible(gc())
  }

  ##########################
  #        Figure 3C       #
  ##########################

  fig_name <- paste0("figure-3C", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Suitability Units Figure.", color = "indicator")
  } else {
    prob_stacks <- c("totalMean", "totalMedian", "totalMax")
    prob_list <- list(length(prob_stacks))

    for (i in 1:length(prob_stacks)) {
      stack_val <- prob_stacks[i]

      prob_vals <- parallel_spec_handler(
        spec.dirs = sp_dirs,
        dir = paste0(log_dir, "/suitability"),
        hv.project.method = "probability",
        col.n = "species-9",
        out.order = paste0("-", stack_val),
        fun = get_prob_stats,
        verbose = verbose
      )

      vals_stack <- stack_projections(
        filenames = prob_vals$filename,
        projection = out_projection,
        projection.method = "bilinear",
        out.dir = paste0(log_dir, "/stack-suitability-", stack_val),
        verbose = verbose
      )

      prob_list[[i]] <- vals_stack
    }

    world_map <- get_world_map(projection = out_projection)

    visualize_suit_units(
      stack.mean = prob_list[[1]],
      stack.median = prob_list[[2]],
      stack.max = prob_list[[3]],
      region = world_map,
      region.name = vis.region.name,
      extent = region_ext,
      projection = out_projection,
      vis.wrap = 1,
      vis.title = vis.title,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      return = TRUE,
      plot.show = plot.show,
      verbose = verbose
    )

    rm(prob_list, world_map)
    invisible(gc())
  }

  ##########################
  #        Figure 3D       #
  ##########################

  fig_name <- paste0("figure-3D-7-9", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Potential Area of Occupancy x Suitability Figure.", color = "indicator")
  } else {
    paoo_file <- paste0(log_dir, "/area-of-occupancy/0.5-inclusion/area-of-occupancy.csv")

    paoo_files <- parallel_spec_handler(
      spec.dirs = sp_dirs,
      shape = shape,
      dir = dirname(dirname(paoo_file)),
      hv.project.method = "0.5-inclusion",
      col.n = "species-9",
      out.order = "-TPAoO",
      fun = get_paoo,
      verbose = verbose
    )

    paoo_files_subset <- copy(paoo_files$filename)

    prob_paoo_files <- gsub("0.5-inclusion", "probability", paoo_files_subset)

    figure_names <- c("1-3", "4-6", "7-9")
    batch <- length(prob_paoo_files) / length(figure_names)

    for (i in 1:batch) {
      start_i <- (i - 1) * batch + 1
      end_i <- i * batch

      prob_batch <- prob_paoo_files[start_i:end_i]
      paoo_batch <- paoo_files_subset[start_i:end_i]
      f_names <- figure_names[i]

      paoo_files <- stack_projections(
        filenames = paoo_batch,
        projection = out_projection,
        projection.method = "near",
        out.dir = paste0(log_dir, "/stack-aoo"),
        binary = TRUE,
        verbose = verbose
      )

      prob_files <- stack_projections(
        filenames = prob_batch,
        projection = out_projection,
        projection.method = "bilinear",
        out.dir = paste0(log_dir, "/stack-prob-aoo"),
        binary = FALSE,
        verbose = verbose
      )

      world_map <- get_world_map(projection = out_projection)

      visualize_dist_suit(
        stack.paoo = paoo_files,
        stack.suitability = prob_files,
        region = world_map,
        region.name = vis.region.name,
        extent = region_ext,
        projection = out_projection,
        vis.wrap = 1,
        vis.title = vis.title,
        save.dir = plot_dir,
        save.name = paste0("figure-3D-", f_names),
        save.device = vis.save.device,
        save.unit = vis.save.unit,
        verbose = verbose
      )
    }

    rm(paoo_files, prob_files, world_map)
    invisible(gc())
  }

  ##########################
  #        Figure 4        #
  ##########################

  fig_name <- paste0("figure-4A", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Composition Figure.", color = "indicator")
  } else {
    paoo_file <- paste0(result_dir, "/area-of-occupancy.csv")
    # Figure 4: stacked barplot wtih taxa per sub region

    richness_dt <- get_taxon_richness(
      paoo.file <- paoo_file,
      stats.file <- sp_stats,
      vis.composition.taxon
    )

    visualize_composition(
      dt = richness_dt,
      region.name = vis.region.name,
      vis.x = "subRegionName",
      vis.x.sort = "westEast",
      vis.y = "relativeRichness",
      vis.fill = vis.composition.taxon,
      vis.group = "group",
      vis.title = vis.title,
      save.dir = plot_dir,
      save.device = vis.save.device,
      verbose = verbose
    )

    richness_write <- richness_dt[, .(get(vis.composition.taxon), relativeRichness, subRegionName, country, groupRelativeRichness, group)]

    setnames(richness_write, "V1", vis.composition.taxon)

    fwrite(richness_write, paste0(result_dir, "/taxon-composition.csv"), bom = TRUE)

    catn("Writing taxonomic composition to markdown file.")

    mdwrite(
      config$files$post_seq_md,
      text = "2;Potential Taxonomic Composition"
    )

    top_orders_md <- richness_dt[, .(get(vis.composition.taxon), relativeRichness, subRegionName, group, groupRelativeRichness, westEast)]

    top_orders_md <- top_orders_md[order(-groupRelativeRichness), head(.SD, 3), by = .(group, subRegionName)]

    setorder(top_orders_md, subRegionName)
    setorder(top_orders_md, westEast)

    top_orders_md <- top_orders_md[, .(V1, relativeRichness, subRegionName, group)]

    setnames(top_orders_md, "V1", vis.composition.taxon)

    mdwrite(
      config$files$post_seq_md,
      text = "3;Top 3 per Group in Each Sub-Region",
      data = top_orders_md
    )

    richness_group_md <- richness_dt[, .(group, groupRelativeRichness, subRegionName, westEast)]
    setorder(richness_group_md, subRegionName)
    setorder(richness_group_md, westEast)
    richness_group_md <- richness_group_md[, westEast := NULL]
    richness_group_md <- unique(richness_group_md, by = c("group", "subRegionName"))

    mdwrite(
      config$files$post_seq_md,
      text = "3;Taxonomic Groups",
      data = richness_group_md
    )

    rm(paoo_file, richness_dt, top_orders_md, richness_group_md)
    invisible(gc())
  }

  ##########################
  #        Figure 5        #
  ##########################

  fig_name <- paste0("figure-5C", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Connections figure.", color = "indicator")
  } else {
    paoo_file <- paste0(result_dir, "/area-of-occupancy.csv")
    # Figure 4: stacked barplot wtih taxa per sub region

    richness_dt <- get_taxon_richness(
      paoo.file <- paoo_file,
      stats.file <- sp_stats,
      "species"
    )

    taxon_names <- c("species", "family", "order")
    figure_names <- c("A", "B", "C")
    connections_md <- copy(richness_dt)

    for (i in 1:length(taxon_names)) {
      taxon <- taxon_names[i]
      figure <- figure_names[i]

      connections_dt <- get_connections(
        dt = richness_dt,
        taxon = taxon,
        verbose = verbose
      )

      # Figure 5: Sankey with floristic regions
      visualize_connections(
        dt = connections_dt,
        taxon = taxon,
        region.name = vis.region.name,
        subregion.name = vis.subregion.name,
        vis.title = vis.title,
        save.dir = plot_dir,
        save.device = vis.save.device,
        save.unit = vis.save.unit,
        save.name = paste0("figure-5", figure),
        plot.show = plot.show,
        verbose = verbose
      )
    }

    connections_md <- get_connections(
      dt = richness_dt,
      taxon = "species",
      verbose = verbose
    )

    subsets <- c("subRegionName", "country", "originCountry", "originCountryCode", "connections")

    con_spec <- connections_md[, c(.SD, mget(subsets)), .SDcols = "cleanName"]

    fwrite(con_spec, paste0(result_dir, "/connections-species.csv"), bom = TRUE)

    con_fam <- connections_md[, c(.SD, mget(subsets)), .SDcols = "family"]

    con_fam <- con_fam[, connections := sum(connections, na.rm = TRUE), by = c("family", "subRegionName")]

    con_fam <- unique(con_fam, by = "family")

    fwrite(con_fam, paste0(result_dir, "/connections-family.csv"), bom = TRUE)

    con_order <- connections_md[, c(.SD, mget(subsets)), .SDcols = "order"]

    con_order <- con_order[, connections := sum(connections, na.rm = TRUE), by = c("order", "subRegionName")]

    con_order <- unique(con_order, by = "order")

    fwrite(con_order, paste0(result_dir, "/connections-order.csv"), bom = TRUE)


    con_subregion <- copy(connections_md)

    con_subregion <- con_subregion[, totalConnections := sum(connections, na.rm = TRUE), by = c("family", "subRegionName")]

    con_subregion <- con_subregion[, .(totalConnections, subRegionName)]

    con_subregion <- unique(con_subregion, by = "subRegionName")

    setorder(con_subregion, -totalConnections)

    mdwrite(
      config$files$post_seq_md,
      text = "2;Potential Family Connection Counts in each Sub-Region",
      data = con_subregion
    )

    con_origin <- copy(connections_md)

    con_origin <- con_origin[, totalConnections := sum(connections, na.rm = TRUE), by = c("family", "originCountry")]

    con_origin <- con_origin[, .(totalConnections, originCountry)]

    con_origin <- unique(con_origin, by = "originCountry")

    setorder(con_origin, -totalConnections)

    mdwrite(
      config$files$post_seq_md,
      text = "2;Potential Connection Counts from each Origin Country",
      data = con_origin
    )

    rm(richness_dt, connections_md, con_spec, con_fam, con_order, subsets, con_subregion, con_origin)
    invisible(gc())
  }

  fig_name <- paste0("figure-6-log", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Species Latitudinal Ranges figure.", color = "indicator")
  } else {
    lat_stats <- copy(sp_stats)
    lat_stats <- lat_stats[, .(cleanName, order, overlapRegion, observations)]
    lat_stats <- get_order_group(
      lat_stats
    )
    spec_name_stats <- unique(lat_stats$cleanName)

    spec_list <- readLines(paste0(hv.dir, "/", hv.method, "-sequence/stats/spec-iteration-list.txt"))

    spec_names <- gsub("-", " ", gsub(".csv", "", basename(spec_list)))

    spec_list <- spec_list[spec_names %in% spec_name_stats]

    spec_name_stats <- unique(lat_stats, by = c("cleanName", "overlapRegion"))

    spec_count_dt <- count_observations(
      spec.list = spec_list,
      dimensions = 4,
      method = "median",
      verbose = verbose
    )

    setnames(spec_count_dt, "species", "cleanName")

    spec_count_dt <- spec_count_dt[, .(cleanName, medianLat)]

    catn("Species min absolute median latitude:", highcat(min(spec_count_dt$medianLat)))

    catn("Species max absolute median latitude:", highcat(max(spec_count_dt$medianLat)))

    merged_dt <- merge(lat_stats, spec_count_dt, by = "cleanName")

    merged_dt <- unique(merged_dt, by = "cleanName")

    visualize_lat_distribution(
      input.dt = merged_dt,
      model.scale = "log",
      region.name = vis.region.name,
      vis.title = vis.title,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      plot.show = plot.show,
      verbose = verbose
    )

    lat_dist_md <- copy(merged_dt)

    model <- lm(data = lat_dist_md, overlapRegion ~ medianLat)

    model_summary <- summary(model)

    mdwrite(
      config$files$post_seq_md,
      text = "2;Species Latitudinal Ranges"
    )

    post_dir <- "./outputs/post-process/images"
    create_dir_if(post_dir)

    model_md <- model_to_md(model)

    mdwrite(
      config$files$post_seq_md,
      text = model_md,
    )

    mdwrite(
      config$files$post_seq_md,
      text = "3;Residuals vs Fitted",
      image = model,
      image.which = 1,
      image.out = paste0(post_dir, "/residuals-fitted.jpeg")
    )

    mdwrite(
      config$files$post_seq_md,
      text = "3;Q-Q Residuals",
      image = model,
      image.which = 2,
      image.out = paste0(post_dir, "/qq-residuals.jpeg")
    )

    mdwrite(
      config$files$post_seq_md,
      text = "3;Scale-Location",
      image = model,
      image.which = 3,
      image.out = paste0(post_dir, "/scale-location.jpeg")
    )

    mdwrite(
      config$files$post_seq_md,
      text = "3;Residuals vs Leverage",
      image = model,
      image.which = 4,
      image.out = paste0(post_dir, "/residuals-leverage.jpeg")
    )

    rm(lat_dist_md, merged_dt, spec_count_dt, lat_stats)
    invisible(gc())
  }

  # fig_name <- paste0("figure-7",".", vis.save.device)
  # if (fig_name %in% existing_plots) {
  #   vebcat("Skipping Species Latitudinal Ranges figure.", color = "indicator")
  # } else {
  #
  #   filter_dir <- "./outputs/filter"
  #
  #
  #   exclude_dirs <- c("chunk","test","logs")
  #
  #   exclude_files <- c("occ.csv",".zip","logs","chunk")
  #
  #   filter_files <- get_repository_files("filter")
  #
  #   calc_list_rows(filter_files)
  #
  #   filter_files <- get_repository_files("setup", subset = "wrangle")
  #
  #   calc_list_rows(filter_files, begin = 1, end = 3)
  #
  #   progressive_dirname(filter_files, begin = 1, end = NULL)
  #
  # }
}
