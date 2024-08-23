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
  vis.region.name = "The Arctic"
  vis.subregion.name = "Floristic Province"
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
  log_dir_cell <- paste0(log_dir, "/cell")
  log_dir_aoo <- paste0(log_dir, "/area-of-occupancy")
  log_stacks_aoo <- paste0(log_dir_aoo, "-stacks")
  log_dir_suitability <- paste0(log_dir, "/suitability")
  log_stacks_suitability <- paste0(log_dir_suitability, "-stacks")

  # create dirs if they do not exist
  create_dir_if(stats_dir, log_dir, result_dir, rast_dir, plot_dir, log_dir_cell, log_dir_aoo, log_stacks_aoo, log_dir_suitability)

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
  region_ext <- ext(region)

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
  sp_group_dirs <- split_spec_by_group(sp_dirs, sp_stats, "cleanName", is.file = TRUE, verbose)
  
  group_richness <- list()
  
  for (group in names(sp_group_dirs)) {
    group_dirs <- sp_group_dirs[[group]]$filename
    # Get the cell richness per group
    inc_dt <- parallel_spec_handler(
      spec.dirs = group_dirs,
      dir = paste0(log_dir_cell, "/", group),
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
  
  
  all_richness <- group_richness
  all_richness$all <- summed_richness
  
  rm(summed_richness, group_richness)
  invisible(gc())
  

  ##########################
  #        Figure 1        #
  ##########################

  fig_name <- paste0("figure-1A", ".", vis.save.device)
  if (fig_name %in% existing_plots) {
    vebcat("Skipping Frequency Figures.", color = "indicator")
  } else {
    freq_stack <- copy(all_richness$all)
    freq_stack <- freq_stack[!is.na(freq_stack$cellRichness)]
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
    hotspots <- list()
    
    for (richness in names(all_richness)) {
      
      hotspot_raster <- convert_template_raster(
        input.values = all_richness[[richness]]$cellRichness,
        template.filename = template_file,
        projection = config$projection$out,
        projection.method = "near",
        out.dir = paste0(rast_dir, "/", richness),
        verbose = verbose
      )
      
      hotspots[richness] <- hotspot_raster 
    }
    
    world_map <- get_world_map(projection = config$projection$out)
    
    for (hotspot in names(hotspots)) {
      hotspot_raster <- hotspots[[hotspot]]
      
      visualize_hotspots(
        rast = hotspot_raster,
        region = world_map,
        extent = region_ext,
        region.name = vis.region.name,
        projection = config$projection$out,
        save.dir = plot_dir,
        save.name = paste0("figure-2-", hotspot),
        save.device = vis.save.device,
        save.unit = vis.save.unit,
        vis.title = vis.title,
        verbose = verbose
      )
    }


    res_hotspots_nums <- all_richness$all[, .SD[which.max(cellRichness)], by = subRegionName]
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
    paoo_files <- list()
    
    for (group in names(sp_group_dirs)) {
      group_dirs <- sp_group_dirs[[group]]$filename
      
      paoo <- parallel_spec_handler(
        spec.dirs = group_dirs,
        shape = shape,
        dir = paste0(log_dir_aoo, "/", group),
        hv.project.method = "0.5-inclusion",
        col.n = "species-9",
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_stacks <- list()
    
    for (i in 1:(length(paoo_files) + 1)) {
      
      if (i == length(paoo_files) + 1) {
        group_name <- "all"
        group <- combine_top_groups(paoo_files, "-TPAoO")
      } else {
        group_name <- names(paoo_files)[i]
        group <- paoo_files[[i]]
      }
      
      filenames <- group$filename
      
      paoo <- stack_projections(
        filenames = filenames,
        projection = config$projection$out,
        projection.method = "near",
        out.dir = paste0(log_stacks_aoo, "/", group_name),
        binary = TRUE,
        verbose = verbose
      )  
      
      paoo_stacks[[group_name]] <- paoo
    }
    
    world_map <- get_world_map(projection = config$projection$out)

    for (stack in names(paoo_stacks)) {
      visualize_paoo(
        rast = paoo_stacks[[stack]],
        region = world_map,
        region.name = vis.region.name,
        extent = region_ext,
        projection = config$projection$out,
        vis.wrap = 3,
        save.dir = plot_dir,
        save.name = paste0("figure-3A-", stack),
        save.device = vis.save.device,
        save.unit = vis.save.unit,
        vis.title = vis.title,
        verbose = verbose
      )
    }
    
    paoo_files$all <- combine_top_groups(paoo_files, "-TPAoO")

    paoo_md <- paoo_files$all[, 1:3]

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
      dir = log_dir_suitability,
      hv.project.method = "probability",
      col.n = "species-9",
      out.order = "-totalMean",
      fun = get_prob_stats,
      verbose = verbose
    )
    
    all_probs <- fread(paste0(result_dir, "/suitability/suitability.csv"))
    
    split_prob <- split_spec_by_group(all_probs, sp_stats, "cleanName", verbose)
    
    split_prob <- lapply(split_prob, function(x) {
      x <- unique(x, by = "species")
      setorder(x, -totalMean)
      x <- x[1:9]
      setcolorder(x, c(setdiff(names(x), "filename"), "filename"))
    })
    
    
    prob_stacks <- list()
    
    for (i in 1:(length(split_prob) + 1)) {
      if (i == length(split_prob) + 1) { # +1 for all together to save memory
        group_name <- "all"
        group <- combine_top_groups(split_prob, "-totalMean")
      } else {
        group <- split_prob[[i]]
        group_name <- names(split_prob)[i]
      }
      
      prob_stack <- stack_projections(
        filenames = group$filename,
        projection = config$projection$out,
        projection.method = "bilinear",
        out.dir = paste0(log_stacks_suitability, "/totalMean/", group_name),
        verbose = verbose
      )
      
      prob_stacks[[group_name]] <- prob_stack
    }
    
    world_map <- get_world_map(projection = config$projection$out)
    
    for (stack in names(prob_stacks)) {
      visualize_suitability(
        stack = prob_stacks[[stack]],
        region = world_map,
        region.name = vis.region.name,
        extent = region_ext,
        projection = config$projection$out,
        vis.unit = "mean",
        vis.wrap = 3,
        vis.title = vis.title,
        save.dir = plot_dir,
        save.name = paste0("figure-3B-", stack),
        save.device = vis.save.device,
        save.unit = vis.save.unit,
        verbose = verbose
      )
    }

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
        shape = shape,
        dir = log_dir_suitability,
        hv.project.method = "probability",
        col.n = "species-9",
        out.order = paste0("-", stack_val),
        fun = get_prob_stats,
        verbose = verbose
      )

      vals_stack <- stack_projections(
        filenames = group$filename,
        projection = config$projection$out,
        projection.method = "bilinear",
        out.dir = paste0(log_stacks_suitability, "/", stack_val, "/all"),
        verbose = verbose
      )

      prob_list[[i]] <- vals_stack
    }

    world_map <- get_world_map(projection = config$projection$out)

    visualize_suit_units(
      stack.mean = prob_list[[1]],
      stack.median = prob_list[[2]],
      stack.max = prob_list[[3]],
      region = world_map,
      region.name = vis.region.name,
      extent = region_ext,
      projection = config$projection$out,
      vis.wrap = 1,
      vis.title = vis.title,
      save.dir = plot_dir,
      save.name = "figure-3C",
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
    paoo_files <- list()
    
    for (group in names(sp_group_dirs)) {
      group_dirs <- sp_group_dirs[[group]]$filename
      
      paoo <- parallel_spec_handler(
        spec.dirs = group_dirs,
        shape = shape,
        dir = paste0(log_dir_aoo, "/", group),
        hv.project.method = "0.5-inclusion",
        col.n = "species-9",
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_files_subset <- combine_top_groups(paoo_files, "-TPAoO")

    paoo_files_subset <- copy(paoo_files_subset$filename)

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
        projection = config$projection$out,
        projection.method = "near",
        out.dir = paste0(log_stacks_aoo, "/top-", f_names),
        binary = TRUE,
        verbose = verbose
      )

      prob_files <- stack_projections(
        filenames = prob_batch,
        projection = config$projection$out,
        projection.method = "bilinear",
        out.dir = paste0(log_stacks_suitability, "/totalMean/top-", f_names),
        binary = FALSE,
        verbose = verbose
      )

      world_map <- get_world_map(projection = config$projection$out)

      visualize_dist_suit(
        stack.paoo = paoo_files,
        stack.suitability = prob_files,
        region = world_map,
        region.name = vis.region.name,
        extent = region_ext,
        projection = config$projection$out,
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
    # Figure 4: stacked barplot wtih taxa per sub region
    paoo_files <- list()
    
    for (group in names(sp_group_dirs)) {
      group_dirs <- sp_group_dirs[[group]]$filename
      
      paoo <- parallel_spec_handler(
        spec.dirs = group_dirs,
        shape = shape,
        dir = paste0(log_dir_aoo, "/", group),
        hv.project.method = "0.5-inclusion",
        col.n = "species-9",
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_file <- combine_top_groups(paoo_files, "-TPAoO")

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
    paoo_files <- list()
    
    for (group in names(sp_group_dirs)) {
      group_dirs <- sp_group_dirs[[group]]$filename
      
      paoo <- parallel_spec_handler(
        spec.dirs = group_dirs,
        shape = shape,
        dir = paste0(log_dir_aoo, "/", group),
        hv.project.method = "0.5-inclusion",
        col.n = "species-9",
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_file <- combine_top_groups(paoo_files, "-TPAoO")
    
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

      # Figure 5: connections map
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
        verbose = TRUE
      )
    }

    connections_md <- get_connections(
      dt = richness_dt,
      taxon = "species",
      verbose = verbose
    )

    subsets <- c("subRegionName", "country", "originCountry", "originCountryCode", "connections")

    con_spec <- connections_md[, c(.SD, mget(subsets)), .SDcols = "cleanName"]
    con_res_dir <- paste0(result_dir, "/connections")
    create_dir_if(con_res_dir)

    fwrite(con_spec, paste0(con_res_dir, "/species.csv"), bom = TRUE)

    con_fam <- connections_md[, c(.SD, mget(subsets)), .SDcols = "family"]

    con_fam <- con_fam[, connections := sum(connections, na.rm = TRUE), by = c("family", "subRegionName")]

    con_fam <- unique(con_fam, by = "family")

    fwrite(con_fam, paste0(con_res_dir, "/family.csv"), bom = TRUE)

    con_order <- connections_md[, c(.SD, mget(subsets)), .SDcols = "order"]

    con_order <- con_order[, connections := sum(connections, na.rm = TRUE), by = c("order", "subRegionName")]

    con_order <- unique(con_order, by = "order")

    fwrite(con_order, paste0(con_res_dir, "/order.csv"), bom = TRUE)


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
    
    merged_dt <- spec_count_dt[lat_stats, on = "cleanName"]

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
