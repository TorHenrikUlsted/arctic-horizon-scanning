source_all("./src/visualize/components")

visualize_sequence <- function(out.dir = "./outputs/visualize", shape, projection = "longlat", hv.dir, hv.method = "box", verbose = FALSE) {
  
  # Test params
  shape = "./outputs/setup/region/cavm-noice/cavm-noice.shp"
  hv.dir = paste0("./outputs/hypervolume/test-small")

  ##########################
  #       Load dirs        #
  ##########################
  log_dir <- paste0(out.dir, "/logs")
  get_vis_log <- paste0(log_dir, "/get-visualize-data")
  hv_proj_dir <- paste0(hv.dir, "/", hv.method, "-sequence/projections")
  
  # create dirs
  create_dir_if(c(log_dir, get_vis_log))
  
  ##########################
  #       Load files       #
  ##########################
  # Set up error handling
  warn_file <- paste0(log_dir, "/", hv.method, "-warning.txt")
  err_file <- paste0(log_dir, "/", hv.method, "-error.txt")
  stats_file <- paste0(hv.dir, "/", hv.method, "-sequence/stats.csv")
  
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
  sp_stats <- get_stats_data(hv.dir, hv.method, verbose = F, warn = warn, err = err)
  
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
    verbose = verbose
  )
  
  # Get number of cells per floristic region to calculate proportional values to each region and make them comparable
  freq_stack <- na.omit(inc_dt)
  
  freq_stack <- freq_stack[, ncellsRegion := .N, by = floregName]
  
  freq_stack <- freq_stack[, propCells := ncellsRegion / nrow(freq_stack), by = floregName]
  
  visualize_freqpoly(
    sp_cells = freq_stack, 
    region = region, 
    region.name = "Arctic", 
    vis.x = "richness",
    vis.title = TRUE,
    vis.color = "country", 
    vis.shade = "floregName",
    vis.shade.name = "FloristicProvince",
    verbose = verbose
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
    projection = out_projection
  )
  
  world_map <- get_world_map(projection = out_projection)
 
  visualize_hotspots(
    rast = hotspot_raster,
    region = world_map,
    extent = region_ext,
    region.name = "Arctic",
    projection  = out_projection,
    projection.method = "near",
    vis.title = TRUE
  )
  
  inc_cover <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/coverage"), 
    n = 9,
    out.order = "coverage",
    fun = get_inc_coverage,
    verbose = FALSE
  )
  
  inc_cover <- stack_projections(
    filenames = inc_cover, 
    projection = out_projection, 
    projection.method = "near", 
    binary = TRUE
  )
  
  visualize_highest_spread(
    rast = inc_cover,
    region = world_map,
    region.name = "Arctic",
    extent = region_ext,
    projection  = out_projection,
    projection.method = "near",
    verbose = FALSE
  )
  
  ##########################
  #        Figure 3        #
  ##########################
  
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  # Get the max values of all probability raster files
  prob_max <- parallel_spec_handler(
    spec.dirs = sp_dirs,
    dir = paste0(get_vis_log, "/max"),
    hv.project.method = "probability",
    n = 3,
    out.order = "maxValue",
    fun = get_prob_max,
    verbose = FALSE
  )
  
  # Read the rasters with the highest max values
  prob_stack <- stack_projections(
    filenames = prob_max, 
    projection = out_projection, 
    projection.method = "bilinear", 
  )
  
  visualize_suitability(
    rast = prob_stack, 
    region = region, 
    region.name = "Arctic"
  )
  
  ##########################
  #        Figure 4        #
  ##########################
  # Figure 4: stacked barplot wtih taxa per floristic region
  region_richness_dt <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/region-richness"),
    shape = shape,
    hv.project.method = "inclusion-0.5", 
    fun = get_region_richness, 
    verbose = FALSE
  )
  
  # Append taxaon rank
  taxon_richness <- append_taxon(region_richness_dt, "species", verbose = TRUE)
  
  # Calculate richness for specified taxon
  richness_dt <- calculate_taxon_richness(taxon_richness, "order")
  
  richness_dt <- handle_region_dt(richness_dt)
  
  print(head(richness_dt, 3))
  
  visualize_richness(
    dt = richness_dt,
    region.name = "Arctic",
    vis.x = "floristicProvince", 
    vis.x.sort = "sn",
    vis.y = "relativeRichness", 
    vis.fill = "order",
    vis.group = "order",
    verbose = FALSE
  )
  
  
  ##########################
  #        Figure 5        #
  ##########################
  
  sankey_dt <- append_src_country(sp_stats, richness_dt, level = "lvl2Name", taxon = "species", verbose = FALSE)
  # Figure 5: Sankey with floristic regions 
  visualize_sankey(
    dt = sankey_dt,
    taxon = "species",
    level = "lvl2Name",
    verbose = FALSE
  )

  # Figure 5: Sankey with floristic regions 
  visualize_connections(
    dt = sankey_dt,
    verbose = FALSE
  )
  
}
