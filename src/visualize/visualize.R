source_all("./src/visualize/components")

visualize_sequence <- function(spec.list, out.dir, hv.dir, hv.method, x.threshold, verbose) {
  
  # Test params
  spec.list = list.files("./outputs/filter/arctic/chunk/species", full.names = TRUE)
  out.dir = "./outputs/visualize" 
  hv.dir = "./outputs/hypervolume/sequence"
  hv.method = "box"
  hv.projection = longlat_crs
  x.threshold = 0.2
  show.plot = FALSE
  verbose = T

  ##########################
  #       Load dirs        #
  ##########################
  log_dir <- paste0(out.dir, "/logs")
  get_vis_log <- paste0(log_dir, "/get-visualize-data")
  
  # create dirs
  create_dir_if(c(log_dir, get_vis_log))
  
  ##########################
  #       Load files       #
  ##########################
  # Set up error handling
  warn_file <- paste0(log_dir, "/", hv.method, "-warning.txt")
  err_file <- paste0(log_dir, "/", hv.method, "-error.txt")
  
  create_file_if(c(warn_file, err_file))
  
  ##########################
  #     Load varaibles     #
  ##########################
  
  sp_dirs <- list.dirs(paste0(hv.dir, "/projections/", hv.method))
  
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]
  
  cavm <- load_region("./resources/region/cavm2003/cavm.shp")
  cavm <- load_region("./outputs/setup/region/cavm-noice/cavm-noice.shp")
  
  cavm <- handle_region(cavm)
  
  cavm_laea <- terra::project(cavm, laea_crs)
  
  cavm_laea_ext <- ext(cavm_laea)
  
  ##########################
  #       Check data       #
  ##########################
  source_all("./src/visualize/components")
  # Check the output if it is in the expected format and length
  checked_data <- check_hv_output(spec.list, hv.dir, hv.method, hv.projection, hv.inc.t = 0.5, verbose = FALSE)
  
  # Clean the output data - 9 hours
  cleaned_data <- clean_hv_output(checked_data, out.dir, hv.dir, hv.projection, verbose = T)
  
  # Get stats csv file
  sp_stats <- get_stats_data(cleaned_data, hv.dir, hv.method, verbose = F, warn = warn, err = err)
  
  ##########################
  #        Figure 1        #
  ##########################
  
  inc_dt <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/cell"), 
    region = shapefiles[[1]],
    hv.project.method  = "inclusion-0.5",
    test = 0,
    fun = get_inclusion_cell,
    batch = TRUE,
    node.log.append = FALSE,
    verbose = FALSE
  )
  
  
  # Get number of cells per floristic region to calculate proportional values to each region and make them comparable
  freq_stack <- na.omit(inc_dt)
  
  freq_stack <- freq_stack[, ncellsRegion := .N, by = floristicProvince]
  
  freq_stack <- freq_stack[, propCells := ncellsRegion / nrow(freq_stack), by = floristicProvince]
  
  visualize_freqpoly(
    sp_cells = freq_stack, 
    region = cavm, 
    region.name = "Arctic", 
    vis.x = "richness",
    vis.title = TRUE,
    vis.color = "country", 
    vis.shade = "floristicProvince",
    verbose = FALSE
  )
  
  ##########################
  #        Figure 2        #
  ##########################
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  
  # Convert to raster by using a template
  hotspot_raster <- convert_template_raster(
    input.values = inc_dt$richness,
    hv.dir = hv.dir, 
    hv.method = hv.method,
    hv.project.method = "inclusion-0.5",
    projection = laea_crs
  )
  
  world_map <- get_world_map(projection = laea_crs)
  
 
  visualize_hotspots(
    rast = hotspot_raster,
    region = world_map,
    extent = cavm_laea_ext,
    region.name = "Arctic",
    projection  = laea_crs,
    projection.method = "near",
    vis.title = TRUE
  )
  
  inc_cover <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/coverage"), 
    region = NULL,
    n = 9,
    out.order = "coverage",
    fun = get_inc_coverage,
    verbose = FALSE
  )
  
  inc_cover <- stack_projections(
    filenames = inc_cover, 
    projection = laea_crs, 
    projection.method = "near", 
    binary = TRUE
  )
  source_all("./src/visualize/components")
  visualize_highest_spread(
    rast = inc_cover,
    region = world_map,
    region.name = "Arctic",
    extent = cavm_laea_ext,
    projection  = laea_crs,
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
    region = NULL,
    hv.project.method  = "probability",
    n = 3,
    out.order = "maxValue",
    fun = get_prob_max,
    verbose = FALSE
  )
  
  # Read the rasters with the highest max values
  prob_stack <- stack_projections(
    filenames = prob_max, 
    projection = laea_crs, 
    projection.method = "bilinear", 
  )
  
  visualize_suitability(
    rast = prob_stack, 
    region = cavm_laea, 
    region.name = "Arctic"
  )
  
  ##########################
  #        Figure 4        #
  ##########################
  # Figure 4: stacked barplot wtih taxa per floristic region
  region_richness_dt <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/region-richness"),
    region = shapefiles[[1]],
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
  source_all("./src/visualize/components")
  # Figure 5: Sankey with floristic regions 
  visualize_sankey(
    dt = sankey_dt,
    taxon = "species",
    level = "lvl2Name",
    verbose = FALSE
  )
  
}
