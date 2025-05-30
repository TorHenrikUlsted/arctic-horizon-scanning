source_all("./src/visualize/components")

visualize_sequence <- function(out.dir = "./outputs/visualize", res.unknown, res.known, shape, hv.dir, hv.method = "box", vis.projection = "longlat", vis.title = TRUE, vis.region.name = "Region", vis.subregion.name = "Sub Region", vis.composition.taxon = "order", vis.save.device = "jpeg", vis.save.unit = "px", validation = FALSE, plot.show = FALSE, verbose = FALSE) {
  
  #------------------------#
  ####    Load dirs     ####
  #------------------------#
  
  hv_dir <- file.path(hv.dir,  paste0(hv.method, "-sequence"))
  hv_proj_dir <- file.path(hv_dir, "projections")
  hv_stats_dir <- file.path(hv_dir, "stats")
  vis_dir <- file.path(out.dir, paste0(hv.method, "-sequence"))
  stats_dir <- file.path(vis_dir, "stats")
  log_dir <- file.path(vis_dir, "logs")
  result_dir <- file.path(vis_dir, "results")
  rast_dir <- file.path(log_dir, "rasters")
  plot_dir <- file.path(vis_dir, "plots")
  log_dir_cell <- file.path(log_dir, "cell")
  log_dir_aoo <- file.path(log_dir, "area-of-occupancy")
  log_stacks_aoo <- paste0(log_dir_aoo, "-stacks")
  log_dir_suitability <- file.path(log_dir, "suitability")
  log_stacks_suitability <- paste0(log_dir_suitability, "-stacks")
  
  # create dirs if they do not exist
  create_dir_if(stats_dir, log_dir, result_dir, rast_dir, plot_dir, log_dir_cell, log_dir_aoo, log_stacks_aoo, log_dir_suitability)
  
  #------------------------#
  ####   Load files     ####
  #------------------------#
  # Set up error handling
  warn_file <- file.path(log_dir, paste0(hv.method, "-warning.txt"))
  err_file <- file.path(log_dir, paste0(hv.method, "-error.txt"))
  stats_file <- file.path(hv.dir, paste0(hv.method, "-sequence/stats/stats.csv"))
  hv_removed_sp <- file.path(hv.dir, paste0(hv.method, "-sequence/stats/removed-species.csv"))
  res_file <- file.path(hv.dir, paste0(hv.method, "-sequence/stats/stats.csv"))
  
  expected_res_file <- file.path("./outputs/filter", gsub("_", "-", res.unknown), "absent-final.csv")
  unknown_chunk_dir <- file.path("./outputs/filter", gsub("_", "-", res.unknown), "chunk/species")
  known_chunk_dir <- file.path("./outputs/filter", gsub("_", "-", res.known), "chunk/species")
  
  known_res_file <- file.path("./outputs/filter", gsub("_", "-", res.known), "present-final.csv")
  
  if (vis.title) {
    plot_dir <- file.path(plot_dir, "title", vis.save.device)
  } else {
    plot_dir <- file.path(plot_dir, "no-title", vis.save.device)
  }
  
  create_file_if(warn_file, err_file)
  
  frequency_files <- c("figure-1A-descriptive", "figure-1A", "figure-1B")
  hotspot_files <- paste0("figure-2-", c("all", "angiosperms", "gymnosperms", "pteridophytes", "combined"))
  occpancy_files <- paste0("figure-3A-", c("all", "angiosperms", "gymnosperms", "pteridophytes"))
  suitability_files <- paste0("figure-3B-", c("all", "angiosperms", "gymnosperms", "pteridophytes"))
  suit_occupancy_files <- "figure-3D-1-5"
  composition_files  <- c("figure-4A", "figure-AB", "figure-4B")
  connections_files <- "figure-5-all"
  gamlss_files <- c("figure-6A", "figure-6B", "figure-6C")
  gam_files <- "figure-7"
  polynomial_files <- "figure-8"
  sankey_files <- "figure-9"
  
  
  # Initial full check
  all_files <- c(
    frequency_files,
    hotspot_files,
    occpancy_files,
    suitability_files,
    suit_occupancy_files,
    composition_files,
    connections_files,
    gamlss_files, 
    gam_files,
    polynomial_files,
    sankey_files
  )
  
  check_result <- check_visualization_files(all_files, plot_dir, vis.title, vis.save.device)
  
  if (check_result$complete) {
    catn("All visualization files exist in:", colcat(plot_dir, color = "output"))
    return(invisible())
  } else {
    catn("Directory:", colcat(plot_dir, color = "output"))
    if (length(check_result$existing) > 0) {
      vebcat("Existing files:", highcat(paste(basename(check_result$existing), collapse = ", ")), veb = verbose)
    }
    catn("Missing files:", highcat(paste(basename(check_result$missing), collapse = ", ")))
    catn("Continuing with visualization sequence for missing files...")
  }
  
  #------------------------#
  ####  Load variables  ####
  #------------------------#
  
  sp_dirs <- list.dirs(hv_proj_dir, full.names = TRUE)
  
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]
  
  region <- load_region(shape)
  region <- check_crs(region, vis.projection)
  region_ext <- ext(region)
  
  #------------------------#
  ####    Check Data    ####
  #------------------------#
  
  # Clean species that are not supposed to be there
  stats_clean <- check_output_species(
    out.dir = stats_dir,
    dt.result = res_file,
    dt.expected = expected_res_file,
    hv.removed = hv_removed_sp,
    clean.dirs = hv_proj_dir,
    verbose = verbose
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
  
  # Calculate the number of occurrence in the region

  occ_files <- file.path(
    if (!validation) unknown_chunk_dir else known_chunk_dir,
    paste0(gsub(" ", config$species$file_separator, unique(sp_stats$included$cleanName)), ".csv")
  )
  
  if (!file.exists(file.path(stats_dir, "region-occ.csv"))) {
    overlap_dt <- loop_occ_overlap(
      spec.occ.dir = occ_files,
      region = shape,
      region.name = "Region",
      file.out = file.path(stats_dir, "region-occ.csv")
    )
    
    overlap_dt <- overlap_dt[, .(species, inRegion)]
    print(overlap_dt)
    
    setkey(overlap_dt, species)
    setkey(sp_stats$included, species)
    
    sp_stats$included[overlap_dt, inRegion := i.inRegion]
    
    fwrite(sp_stats$included, file.path(stats_dir, "included-species.csv"))
  }
  
  catn("Found", highcat(length(sp_dirs)), "projections and", highcat(length(unique(sp_stats$included$cleanName))), "Species.")
  
  if (length(sp_dirs) > length(unique(sp_stats$included$cleanName))) {
    stop("Found more projections than species.")
  } else if (length(sp_dirs) < length(unique(sp_stats$included$cleanName))) {
    stop("Found more species than projections.")
  }
  
  #------------------------#
  ####  Get cell data   ####
  #------------------------#
  
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
  sp_group_dirs <- split_spec_by_group(
    sp_dirs, 
    sp_stats$included, 
    "cleanName", 
    is.file = TRUE, 
    verbose
  )
  
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
  
  catn("Getting summed richness")
  summed_richness <- copy(group_richness[[1]])
  
  for (i in 2:length(group_richness)) {
    group_data <- group_richness[[i]]
    
    summed_richness[group_data, on = "cell", cellRichness :=
                      fifelse(
                        is.na(x.cellRichness) & is.na(i.cellRichness), NA_real_,
                        fifelse(
                          is.na(x.cellRichness), i.cellRichness,
                          fifelse(
                            is.na(i.cellRichness), x.cellRichness,
                            x.cellRichness + i.cellRichness
                          )
                        )
                      )]
  }
  
  
  all_richness <- group_richness
  all_richness$all <- summed_richness
  
  rm(summed_richness, group_richness)
  invisible(gc())
  
  #------------------------#
  ####    Frequency     ####
  #------------------------#
  
  finished <- length(check_visualization_files(frequency_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
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

  #------------------------#
  ####    Hotspots      ####
  #------------------------#
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  finished <- length(check_visualization_files(hotspot_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
    vebcat("Skipping Hotspots Figure.", color = "indicator")
  } else {
    # Convert to raster by using a template
    hotspots <- list()
    
    for (richness in names(all_richness)) {
      hotspot_raster <- convert_template_raster(
        input.values = all_richness[[richness]]$cellRichness,
        template.filename = template_file,
        projection = config$simulation$projection,
        projection.method = "near",
        out.dir = paste0(rast_dir, "/", richness),
        verbose = verbose
      )
      
      hotspots[richness] <- hotspot_raster
    }
    
    world_map <- get_world_map(projection = config$simulation$projection)
    
    hotplot <- list()
    
    for (hotspot in names(hotspots)) {
      hotspot_raster <- hotspots[[hotspot]]
      
      all_f <- !grepl("all", hotspot)
      
      fig <- visualize_hotspots(
        rast = hotspot_raster,
        region = world_map,
        extent = region_ext,
        region.name = vis.region.name,
        projection = config$simulation$projection,
        save.dir = plot_dir,
        save.name = paste0("figure-2-", hotspot),
        save.device = vis.save.device,
        save.unit = vis.save.unit,
        vis.title = if (all_f) TRUE else vis.title,
        verbose = verbose
      )
      
      if (all_f) {
        hotplot[[hotspot]] <- fig
      }
    }
    
    catn("Combining group plots")
    combined_fig <- grid.arrange(
      grobs = hotplot,
      ncol = 1, 
      nrow = 3,
      heights = c(1, 1, 1)
    )
    
    save_ggplot(
      save.plot = combined_fig,
      save.name = "figure-2-combined",
      save.width = 3050,
      save.height = 3000*3,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      vis.title = vis.title,
      verbose = verbose
    )
    
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

  #-------------------------#
  #### Area of occupancy ####
  #-------------------------#
  
  finished <- length(check_visualization_files(occpancy_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
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
        group <- combine_groups(paoo_files, "-TPAoO", 9)
      } else {
        group_name <- names(paoo_files)[i]
        group <- paoo_files[[i]]
      }
      
      filenames <- group$filename
      
      paoo <- stack_projections(
        filenames = filenames,
        projection = config$simulation$projection,
        projection.method = "near",
        out.dir = paste0(log_stacks_aoo, "/", group_name),
        binary = TRUE,
        verbose = verbose
      )
      
      paoo_stacks[[group_name]] <- paoo
    }
    
    world_map <- get_world_map(projection = config$simulation$projection)
    
    for (stack in names(paoo_stacks)) {
      visualize_paoo(
        rast = paoo_stacks[[stack]],
        region = world_map,
        region.name = vis.region.name,
        extent = region_ext,
        projection = config$simulation$projection,
        vis.row = 3,
        vis.col = 3,
        save.dir = plot_dir,
        save.name = paste0("figure-3A-", stack),
        save.device = vis.save.device,
        save.unit = vis.save.unit,
        vis.title = vis.title,
        verbose = verbose
      )
    }
    
    paoo_files$all <- combine_groups(paoo_files, "-TPAoO", 9)
    
    paoo_md <- paoo_files$all[, 1:3]
    
    mdwrite(
      config$files$post_seq_md,
      text = "2;Highest Total Potential Area of Occupancy",
      data = paoo_md,
    )
    
    
    rm(paoo_files, paoo, paoo_md, world_map)
    invisible(gc())
  }
  
  #------------------------#
  ####   Suitability    ####
  #------------------------#
  
  finished <- length(check_visualization_files(suitability_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
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
    
    split_prob <- split_spec_by_group(all_probs, sp_stats$included, "cleanName", verbose)
    
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
        group <- combine_groups(split_prob, "-totalMean", 9)
      } else {
        group <- split_prob[[i]]
        group_name <- names(split_prob)[i]
      }
      
      prob_stack <- stack_projections(
        filenames = group$filename,
        projection = config$simulation$projection,
        projection.method = "bilinear",
        out.dir = paste0(log_stacks_suitability, "/totalMean/", group_name),
        verbose = verbose
      )
      
      prob_stacks[[group_name]] <- prob_stack
    }
    
    world_map <- get_world_map(projection = config$simulation$projection)
    
    for (stack in names(prob_stacks)) {
      visualize_suitability(
        stack = prob_stacks[[stack]],
        region = world_map,
        region.name = vis.region.name,
        extent = region_ext,
        projection = config$simulation$projection,
        vis.unit = "mean",
        vis.row = 3,
        vis.col = 3,
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

  #------------------------#
  ####   Suit + paoo    ####
  #------------------------#
  
  finished <- length(check_visualization_files(suit_occupancy_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
    vebcat("Skipping Potential Area of Occupancy x Suitability Figure.", color = "indicator")
  } else {
    rows <- 5
    
    paoo_files <- list()
    
    for (group in names(sp_group_dirs)) {
      group_dirs <- sp_group_dirs[[group]]$filename
      
      paoo <- parallel_spec_handler(
        spec.dirs = group_dirs,
        shape = shape,
        dir = paste0(log_dir_aoo, "/", group),
        hv.project.method = "0.5-inclusion",
        col.n = paste0("species-", rows),
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_files_subset <- combine_groups(paoo_files, "-TPAoO", rows)
    
    paoo_files_subset <- copy(paoo_files_subset$filename)
    
    prob_paoo_files <- gsub("0.5-inclusion", "probability", paoo_files_subset)
    
    fig_name <- paste0("1-", rows)
  
    prob_batch <- prob_paoo_files[1:rows]
    paoo_batch <- paoo_files_subset[1:rows]
      
    paoo_files <- stack_projections(
      filenames = paoo_batch,
      projection = config$simulation$projection,
      projection.method = "near",
      out.dir = paste0(log_stacks_aoo, "/top-", fig_name),
      binary = TRUE,
      verbose = verbose
    )
    
    prob_files <- stack_projections(
      filenames = prob_batch,
      projection = config$simulation$projection,
      projection.method = "bilinear",
      out.dir = paste0(log_stacks_suitability, "/totalMean/top-", fig_name),
      binary = FALSE,
      verbose = verbose
    )
    
    world_map <- get_world_map(projection = config$simulation$projection)
    
    visualize_dist_suit(
      stack.paoo = paoo_files,
      stack.suitability = prob_files,
      region = world_map,
      region.name = vis.region.name,
      extent = region_ext,
      projection = config$simulation$projection,
      vis.row = rows,
      vis.col = 1,
      vis.title = vis.title,
      save.dir = plot_dir,
      save.name = paste0("figure-3D-", fig_name),
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      verbose = verbose
    )
    
    prob_mean <- parallel_spec_handler(
      spec.dirs = prob_batch,
      shape = shape,
      dir = log_dir_suitability,
      hv.project.method = "probability",
      col.n = paste0("species-", rows),
      out.order = "-totalMean",
      fun = get_prob_stats,
      verbose = verbose
    )
    
    prob_vals_md <- copy(prob_mean)
    prob_vals_md <- prob_vals_md[, filename := NULL]
    
    prob_vals_md <- prob_vals_md[, 1:rows]
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0(
        "2;Highest Potential Climatic Suitability\n\n",
        "See figure 3D for comparison between **potential area of occupancy** and **potential climatic suitability**"
      ),
      data = prob_vals_md
    )
    
    rm(paoo_files, prob_files, prob_vals_md, world_map)
    invisible(gc())
  }

  #------------------------#
  ####   Composition    ####
  #------------------------#
  
  finished <- length(check_visualization_files(composition_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
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
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_file <- combine_groups(paoo_files, "-TPAoO")
    
    richness_dt <- get_taxon_richness(
      paoo.file = paoo_file,
      stats = sp_stats$included,
      taxon = vis.composition.taxon,
      verbose = verbose
    )
    
    unknown_dt <- get_unknown_composition(
      unknown.path = expected_res_file,
      stats = sp_stats$included,
      out.dir = result_dir,
      verbose = verbose
    )
    
    visualize_composition(
      dt = richness_dt,
      dt.comparison = unknown_dt,
      region.name = vis.region.name,
      vis.x = "subRegionName",
      vis.x.sort = "westEast",
      vis.y = "relativeRichness",
      vis.fill = "order",
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
    
    richness_dt <- get_taxon_richness(
      paoo.file = paoo_file,
      stats = sp_stats$included,
      taxon = "species",
      verbose = verbose
    )
    
    fwrite(richness_dt, paste0(result_dir, "/species-taxon-composition.csv"), bom = TRUE)
    
    rm(paoo_file, richness_dt, top_orders_md, richness_group_md)
    invisible(gc())
  }

  #------------------------#
  ####   Connections    ####
  #------------------------#
  
  finished <- length(check_visualization_files(connections_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
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
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_file <- combine_groups(paoo_files, "-TPAoO")
    
    richness_dt <- get_taxon_richness(
      paoo.file = paoo_file,
      stats = sp_stats$included,
      taxon = "species",
      country = TRUE,
      verbose = verbose
    )
    
    richness_dt <- richness_dt[PAoO > 0]
    
    connections_dt <- get_connections(
      dt = richness_dt,
      verbose = verbose
    )
    
    connections_md <- connections_dt
    
    visualize_spec_ranges_by_subregion(
      dt = connections_dt,
      taxon = "species",
      centroid = FALSE,
      multiple = TRUE,
      region = shape,
      region.name = vis.region.name,
      subregion.name = vis.subregion.name,
      projection = "aeqd",
      vis.title = vis.title,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      save.name = "figure-5",
      plot.show = plot.show,
      verbose = FALSE
    )
    
    subsets <- c("subRegionName", "country", "originCountry", "originCountryCode", "connections")
    
    # Get species connections
    con_spec <- connections_md[, c(.SD, mget(subsets)), .SDcols = "cleanName"]
    con_res_dir <- paste0(result_dir, "/connections")
    create_dir_if(con_res_dir)
    
    fwrite(con_spec, paste0(con_res_dir, "/species.csv"), bom = TRUE)
    
    # Get family connections
    con_fam <- connections_md[, c(.SD, mget(subsets)), .SDcols = "family"]
    
    con_fam <- con_fam[, connections := sum(connections, na.rm = TRUE), by = c("family", "subRegionName")]
    
    con_fam <- unique(con_fam, by = "family")
    
    fwrite(con_fam, paste0(con_res_dir, "/family.csv"), bom = TRUE)
    
    # Get order connections
    con_order <- connections_md[, c(.SD, mget(subsets)), .SDcols = "order"]
    
    con_order <- con_order[, connections := sum(connections, na.rm = TRUE), by = c("order", "subRegionName")]
    
    con_order <- unique(con_order, by = "order")
    
    fwrite(con_order, paste0(con_res_dir, "/order.csv"), bom = TRUE)
    
    # Get SubRegion connections
    con_subregion <- copy(connections_md)
    
    con_subregion <- con_subregion[, totalConnections := uniqueN(cleanName), by = "subRegionName"]
    
    con_subregion <- con_subregion[, .(totalConnections, subRegionName)]
    
    con_subregion <- unique(con_subregion, by = "subRegionName")
    
    setorder(con_subregion, -totalConnections)
    
    mdwrite(
      config$files$post_seq_md,
      text = "2;Potential Species Connection Counts in each Sub-Region",
      data = con_subregion
    )
    
    con_origin <- copy(connections_md)
    
    con_origin <- con_origin[, totalConnections := uniqueN(cleanName), by = "originCountry"]
    
    con_origin <- con_origin[, .(totalConnections, originCountry)]
    
    con_origin <- unique(con_origin, by = "originCountry")
    
    setorder(con_origin, -totalConnections)
    
    mdwrite(
      config$files$post_seq_md,
      text = "2;Potential Species Connection Counts from each Origin Country",
      data = con_origin
    )
    
    rm(richness_dt, connections_md, con_spec, con_fam, con_order, subsets, con_subregion, con_origin)
    invisible(gc())
  }

  #--------------------------------#
  ####   Distribution - GAMLSS  ####
  #--------------------------------#
  
  finished <- length(check_visualization_files(gamlss_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished) {
    vebcat("Skipping GAMLSS figure.", color = "indicator")
  } else {
    
    lat_dt <- analyze_latitudinal_pattern(
      file.path(hv_stats_dir, "stats.csv"),
      file.path(hv_stats_dir, "spec-iteration-list.txt"),
      vis_dir
    )
    
    model <- compare_gamlss_models(
      data = lat_dt, 
      response = "overlapRegion", 
      predictor = "medianLat"
    )
    
    print_gamlss_summary(model$models, model$summary)
    
    visualize_gamlss(
      dt = lat_dt,
      model = model,
      region.name = vis.region.name,
      vis.title = vis.title,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      plot.show = plot.show,
      verbose = verbose
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = "2;Species Latitudinal Ranges"
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = "Summary of Selection Criteria for all GAMLSS Models",
      data = model$summary
    )
    
    post_dir <- "./outputs/post-process/images"
    create_dir_if(post_dir)
    
    model_md <- model_to_md(model$models[[1]]) # Always the best model fit
    
    mdwrite(
      config$files$post_seq_md,
      text = model_md,
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = "3;Summary Plot",
      image = model$models[[1]],
      image.out = paste0(post_dir, "/plot-summary.jpeg")
    )
    
    rm(lat_dt, model, model_md)
    invisible(gc())
  }

  #-----------------------------#
  ####   Distribution - GAM  ####
  #-----------------------------#
  
  finished <- length(check_visualization_files(gam_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished)  {
    vebcat("Skipping GAM figure.", color = "indicator")
  } else {
    # Prepare GAM
    
    # Read from hypervolume dir to include excluded species
    lat_dt <- fread(stats_file, select = c("cleanName", "overlapRegion", "level3Name", "level3Lat", "excluded"))
    
    #lat_dt <- sp_stats$included[, .(cleanName, overlapRegion, level3Name, level3Lat)]
    lat_dt <- data.table::setnames(lat_dt, "level3Lat", "centroidLatitude")
    
    print(lat_dt)
    
    model <- compare_gam_models(
      lat_dt, 
      predictor = "centroidLatitude", 
      response = "overlapRegion"
    )
    
    visualize_gam(
      dt = lat_dt, 
      model = model$models, 
      region.name = vis.region.name, 
      vis.gradient = "viridis-b", 
      vis.title = vis.title, 
      save.dir = plot_dir, 
      save.name = "figure-7",
      save.device = vis.save.device, 
      save.unit = vis.save.unit, 
      plot.save = TRUE, 
      plot.show = plot.show, 
      verbose = verbose
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = "2;Species Botanical Centroid Ranges"
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = "Summary of Selection Criteria for all GAM Models",
      data = model$summary
    )
    
    post_dir <- "./outputs/post-process/images"
    create_dir_if(post_dir)
    
    model_md <- model_to_md(model$models[[1]]) # Always the best model fit
    
    mdwrite(
      config$files$post_seq_md,
      text = model_md,
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = "3;Summary Plot",
      image = model$models[[1]],
      image.out = paste0(post_dir, "/plot-summary-GAM.jpeg")
    )
    
    rm(lat_dt, model, model_md)
    invisible(gc())
  }
  
  #------------------------------------#
  ####   Distribution - Polynomial  ####
  #------------------------------------#
  
  finished <- length(check_visualization_files(polynomial_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
  if (finished)  {
    vebcat("Skipping Polynomial quantile regression figure.", color = "indicator")
  } else {
    
    paoo_dt <- list()
    
    for (group in names(sp_group_dirs)) {
      group_dirs <- sp_group_dirs[[group]]$filename
      
      paoo <- parallel_spec_handler(
        spec.dirs = group_dirs,
        shape = shape,
        dir = paste0(log_dir_aoo, "/", group),
        hv.project.method = "0.5-inclusion",
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_dt[[group]] <- paoo
    }
    
    paoo_dt <- combine_groups(paoo_dt, "-TPAoO")
    setnames(paoo_dt, "species", "cleanName")
    paoo_dt <- paoo_dt[, .(cleanName, TPAoO, propTPAoO)]
    # Filter out 0 TPAoO
    #paoo_dt <- paoo_dt[TPAoO > 0]
    paoo_dt <- unique(paoo_dt, by = "cleanName")
    
    print(nrow(paoo_dt))
    
    
    # Calculate centroids and fit model
    analysis_data <- calculate_species_centroids(
      dt = sp_stats$included,
      out.dir = result_dir,
      verbose = verbose
    )
    
    analysis_data[, specMeanCentroidLat := mean(.SD[["speciesCentroidLat"]], na.rm = TRUE), by = "cleanName"]
    
    analysis_data <- unique(analysis_data, by = "cleanName")
    #analysis_data <- analysis_data[excluded == FALSE]
    
    # merge with paoo
    non_model_data <- analysis_data[!(cleanName %in% paoo_dt$cleanName)]
    non_model_data <- non_model_data[, `:=` (
      TPAoO = 0,
      propTPAoO = 0.0000000
    )]
    model_data <- analysis_data[paoo_dt, on = "cleanName", nomatch = NA]
    analysis_data <- rbind(model_data, non_model_data)
    
    analysis_data <- analysis_data[TPAoO > 0]
    
    response <- "overlapRegion"
    predictor <- "specMeanCentroidLat"
    
    analysis_data <- prepare_poly_spec(
      dt = analysis_data,
      response = response, 
      predictor = predictor,
      transform = NULL,
      by.region = FALSE
    )
    
    model_results <- fit_polynomial_model(
      dt = analysis_data,
      response = response, 
      predictor = predictor
    )
    
    visualize_polynomial_model(
      dt = analysis_data,
      model_results = model_results,
      response = response,
      predictor = predictor,
      region.name = "Arctic",
      vis.title = vis.title,
      save.dir = plot_dir,
      plot.show = plot.show,
      verbose = verbose
    )
    
    # Get model summary
    #summary_stats <- summarize_polynomial_model(model_results)
    
    mdwrite(
      source = config$files$post_seq_md,
      text = "2;Polynomial Model Summary",
      data = cli_to_markdown(summarize_polynomial_model(model_results, transform = NULL))
    )
  }
   
  #------------------------#
  ####      Sankey      ####
  #------------------------#
  
  finished <- length(check_visualization_files(sankey_files, plot_dir, vis.title, vis.save.device, verbose)$missing) == 0
  
    if (finished) {
    vebcat("Skipping Species Sankey figure.", color = "indicator")
  } else {
    paoo_files <- list()
    
    for (group in names(sp_group_dirs)) {
      group_dirs <- sp_group_dirs[[group]]$filename
      
      paoo <- parallel_spec_handler(
        spec.dirs = group_dirs,
        shape = shape,
        dir = paste0(log_dir_aoo, "/", group),
        hv.project.method = "0.5-inclusion",
        out.order = "-TPAoO",
        fun = get_paoo,
        verbose = verbose
      )
      
      paoo_files[[group]] <- paoo
    }
    
    paoo_file <- combine_groups(paoo_files, "-TPAoO")
    
    richness_dt <- get_taxon_richness(
      paoo.file = paoo_file,
      stats = sp_stats$included,
      taxon = "species",
      verbose = verbose
    )
    
    richness_dt <- richness_dt[PAoO > 0]
    
    data.table::setnames(richness_dt, c("originCountry", "subRegionName"), c("origin", "destination"))
    
    # Ensure data is complete and sorted
    dt_sank <- richness_dt[!is.na(origin) & !is.na(destination)]
    dt_sank <- dt_sank[order(-relativeRichness)]
    
    # Aggregate the data if needed
    dt_sank <- dt_sank[, .(
      relativeRichness = sum(relativeRichness)
    ), by = c("origin", "destination")]
    
    save_file <- paste0(result_dir, "/sankey.csv")
    
    catn("Saving sankey file to:", colcat(save_file, color = "output"))
    fwrite(dt_sank, save_file, bom = TRUE)
    
    # Figure 8: Sankey map
    visualize_sankey(
      dt = dt_sank,
      taxon = "species",
      region.name = vis.region.name,
      subregion.name = vis.subregion.name,
      vis.title = vis.title,
      save.dir = plot_dir,
      save.device = vis.save.device,
      save.unit = vis.save.unit,
      save.name = "figure-9",
      plot.show = plot.show,
      verbose = Tverbose
    )
  }
  
}
