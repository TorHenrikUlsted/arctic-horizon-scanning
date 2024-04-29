source_all("./src/visualize/components")

visualize_sequence <- function(out.dir = "./outputs/visualize", shape, projection = "longlat", hv.dir, hv.method = "box", vis.title = TRUE, verbose = FALSE) {
  
  # Test params
  out.dir = "./outputs/visualize/glonaf"
  shape = "./outputs/setup/region/cavm-noice/cavm-noice.shp"
  projection = "laea"
  hv.dir = "./outputs/hypervolume/glonaf"
  hv.method = "box"
  verbose = TRUE
  
  ##########################
  #       Load dirs        #
  ##########################
  hv_proj_dir <- paste0(hv.dir, "/", hv.method, "-sequence/projections")
  vis_dir <- paste0(out.dir, "/", hv.method, "-sequence")
  log_dir <- paste0(vis_dir, "/logs")
  get_vis_log <- paste0(log_dir, "/get-visualize-data")
  rast_dir <- paste0(get_vis_log, "/rasters")
  plot_dir <- paste0(vis_dir, "/plots")
  
  # create dirs
  create_dir_if(c(log_dir, rast_dir, plot_dir))
  
  ##########################
  #       Load files       #
  ##########################
  # Set up error handling
  warn_file <- paste0(log_dir, "/", hv.method, "-warning.txt")
  err_file <- paste0(log_dir, "/", hv.method, "-error.txt")
  stats_file <- paste0(hv.dir, "/", hv.method, "-sequence/stats/stats.csv")
  
  create_file_if(c(warn_file, err_file))
  
  ##########################
  #     Load varaibles     #
  ##########################
  
  sp_dirs <- list.dirs(hv_proj_dir, full.names = TRUE)
  
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]
  
  region <- load_region(shape)

  if (projection == "longlat") {
    region <- handle_region(region)
    region_ext <- ext(region)
    out_projection <- longlat_crs
  } else if (projection == "laea") {
    region <- terra::project(region, laea_crs)
    region_ext <- ext(region)
    out_projection <- laea_crs
  }
  
  ##########################
  #       Check data       #
  ##########################
  # Check the output if it is in the expected format and length
  #checked_data <- check_hv_output(sp_dirs, hv.dir, hv.method, hv.projection, hv.inc.t = 0.5, verbose = FALSE)
  
  # Clean the output data - 9 hours
  #cleaned_data <- clean_hv_output(checked_data, out.dir, hv.dir, hv.projection, verbose = T)
  
  # Get stats csv file
  sp_stats <- get_stats_data(vis_dir, hv.dir, hv.method, verbose = verbose, warn = warn, err = err)
  
  catn("Found", highcat(length(sp_dirs)), "species and", highcat(length(unique(sp_stats$cleanName))), "projections.")
  
  if (length(sp_dirs) < length(unique(sp_stats$cleanName))) {
    stop("More projections than species found.")
  } else if (length(sp_dirs) > length(unique(sp_stats$cleanName))) {
    stop("More species than projections found.")
  }
  
  ##########################
  #        Figure 1        #
  ##########################
  
  inc_dt <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/cell"), 
    shape = shape,
    hv.project.method = "0.5-inclusion",
    test = 0,
    fun = get_inclusion_cell,
    batch = TRUE,
    node.log.append = FALSE,
    verbose = FALSE
  )
  
  # Get number of cells per floristic region to calculate proportional values to each region and make them comparable
  freq_stack <- inc_dt[!is.na(inc_dt$richness)]
  
  #freq_stack <- freq_stack[, ncellsRegion := .N, by = floregName]
  
 # freq_stack <- freq_stack[, propCells := ncellsRegion / nrow(freq_stack), by = floregName]
  visualize_freqpoly(
    spec.cells = freq_stack, 
    region = region, 
    region.name = "Arctic", 
    vis.x = "richness",
    vis.title = TRUE,
    vis.color = "country", 
    vis.shade = "floregName",
    vis.shade.name = "FloristicProvince",
    vis.binwidth = 1,
    vis.x.scale = "sqrt",
    vis.peak.threshold = 0.001,
    save.dir = plot_dir,
    save.device = "jpeg",
    verbose = FALSE
  )
  
  ##########################
  #        Figure 2        #
  ##########################
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  
  # Convert to raster by using a template
  hotspot_raster <- convert_template_raster(
    input.values = inc_dt$richness,
    hv.proj.dir = hv_proj_dir, 
    hv.method = hv.method,
    hv.project.method = "0.5-inclusion",
    projection = out_projection,
    out.dir = rast_dir
  )
  
  world_map <- get_world_map(projection = out_projection)
 
  visualize_hotspots(
    rast = hotspot_raster,
    region = world_map,
    extent = region_ext,
    region.name = "Arctic",
    projection  = out_projection,
    vis.gradient = "viridis-B",
    save.dir = plot_dir,
    save.device = "jpeg",
    vis.title = TRUE
  )
  
  ##########################
  #        Figure 3A       #
  ##########################
  
  inc_cover_files <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/coverage"),
    hv.project.method = "0.5-inclusion",
    n = 9,
    out.order = "coverage",
    fun = get_inc_coverage,
    verbose = FALSE
  )
  
  inc_cover <- stack_projections(
    filenames = inc_cover_files$filename, 
    projection = out_projection, 
    projection.method = "near",
    out.dir = paste0(get_vis_log, "/stack-coverage"),
    binary = TRUE,
    verbose = TRUE
  )
  
  dist_fig <- visualize_distributions(
    rast = inc_cover,
    region = world_map,
    region.name = "the Arctic",
    extent = region_ext,
    projection  = out_projection,
    vis.gradient = "viridis-b",
    vis.wrap = 3,
    save.dir = plot_dir,
    save.device = "jpeg",
    vis.title = TRUE,
    verbose = FALSE
  )
  
  ##########################
  #        Figure 3B       #
  ##########################
  
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  # Get the max values of all probability raster files
  prob_vals <- parallel_spec_handler(
    spec.dirs = sp_dirs,
    dir = paste0(get_vis_log, "/suitability"),
    hv.project.method = "probability",
    n = 9,
    out.order = "mean",
    fun = get_prob_stats,
    verbose = verbose
  )
  
  # Read the rasters with the highest max values
  prob_stack <- stack_projections(
    filenames = prob_vals$filename, 
    projection = out_projection, 
    projection.method = "bilinear",
    out.dir = paste0(get_vis_log, "/stack-suitability-", "mean"),
    verbose = verbose
  )
  
  suit_fig <- visualize_suitability(
    stack = prob_stack, 
    region = world_map, 
    region.name = "the Arctic",
    extent = region_ext,
    projection = laea_crs,
    vis.gradient = "viridis-b",
    vis.unit = "mean",
    vis.wrap = 3,
    vis.title = TRUE,
    save.dir = plot_dir,
    save.device = "jpeg",
    verbose = FALSE
  )
  
  
  ##########################
  #        Figure 3C       #
  ##########################
  
  
  prob_median <- parallel_spec_handler(
    spec.dirs = sp_dirs,
    dir = paste0(get_vis_log, "/suitability"),
    hv.project.method = "probability",
    n = 9,
    out.order = "median",
    fun = get_prob_stats,
    verbose = verbose
  )
  
  median_stack <- stack_projections(
    filenames = prob_median$filename, 
    projection = out_projection, 
    projection.method = "bilinear",
    out.dir = paste0(get_vis_log, "/stack-suitability-", "median"),
    verbose = verbose
  )
  
  
  prob_max <- parallel_spec_handler(
    spec.dirs = sp_dirs,
    dir = paste0(get_vis_log, "/suitability"),
    hv.project.method = "probability",
    n = 9,
    out.order = "max",
    fun = get_prob_stats,
    verbose = verbose
  )
  
  max_stack <- stack_projections(
    filenames = prob_max$filename, 
    projection = out_projection, 
    projection.method = "bilinear",
    out.dir = paste0(get_vis_log, "/stack-suitability-", "max"),
    verbose = verbose
  )
  

  visualize_suit_units(
    stack.mean = prob_stack,
    stack.median = median_stack,
    stack.max = max_stack,
    region = world_map, 
    region.name = "the Arctic",
    extent = region_ext,
    projection = laea_crs,
    vis.gradient = "viridis-b",
    vis.wrap = 1,
    vis.title = TRUE,
    save.dir = plot_dir,
    save.device = "jpeg",
    verbose = FALSE
  )
  
  ##########################
  #        Figure 3D       #
  ##########################
  
  prob_cover_files <- gsub("0.5-inclusion", "probability", inc_cover_files$filename)
  
  prob_cover <- stack_projections(
    filenames = prob_cover_files, 
    projection = out_projection, 
    projection.method = "bilinear",
    out.dir = paste0(get_vis_log, "/stack-prob-cover"),
    binary = FALSE,
    verbose = TRUE
  )
  
  visualize_dist_suit(
    stack.distribution = inc_cover,
    stack.suitability = prob_cover,
    region = world_map, 
    region.name = "the Arctic",
    extent = region_ext,
    projection = laea_crs,
    vis.gradient = "viridis-b",
    vis.wrap = 1,
    vis.title = TRUE,
    save.dir = plot_dir,
    save.device = "jpeg",
    verbose = FALSE
  )
  
  
  ##########################
  #        Figure 4        #
  ##########################
  # Figure 4: stacked barplot wtih taxa per floristic region
  
  region_richness_dt <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/region-richness"),
    shape = shape,
    extra = paste0(get_vis_log, "/cell/0.5-inclusion/cell.csv"), # Has to be the same as inc_dt
    hv.project.method = "0.5-inclusion", 
    fun = get_region_richness, 
    verbose = TRUE
  )
  
  setnames(region_richness_dt, "species", "cleanName")
  
  setnames(sp_stats, "country", "originCountry")
  setnames(sp_stats, "countryCode", "originCountryCode")
  setnames(sp_stats, "meanLong", "originMeanLong")
  setnames(sp_stats, "meanLat", "originMeanLat")
  
  merged_dt <- merge(region_richness_dt, sp_stats, by = "cleanName", allow.cartesian = TRUE)
  
  # Calculate richness for specified taxon
  richness_dt <- calculate_taxon_richness(merged_dt, "order")
  
  print(head(richness_dt, 3))
  source_all("./src/visualize/components")
  visualize_richness(
    dt = richness_dt,
    region.name = "the Arctic",
    vis.x = "subRegionName", 
    vis.x.sort = "we",
    vis.y = "relativeRichness", 
    vis.fill = "order",
    vis.group = "order",
    vis.gradient = "viridis-b",
    vis.title = TRUE,
    save.dir = plot_dir,
    save.device = "jpeg",
    verbose = FALSE
  )
  
  
  ##########################
  #        Figure 5        #
  ##########################
  
  connections_dt <- get_connections(
    dt = richness_dt,
    subset = "order",
    verbose = verbose
  )
  
  connections_points <- get_con_points(
    dt = connections_dt,
    verbose = verbose
  )
  
  
  source_all("./src/visualize/components")
  # Figure 5: Sankey with floristic regions 
  visualize_connections(
    dt = connections_dt,
    taxon = "Order",
    region.name = "the Arctic",
    vis.gradient = "viridis-b",
    vis.title = TRUE,
    save.dir = plot_dir,
    save.device = "jpeg",
    verbose = TRUE
  )
  
}
