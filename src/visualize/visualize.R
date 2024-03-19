source_all("./src/visualize/components")

visualize <- function(spec.list, out.dir, hv.dir, hv.method, x.threshold, verbose) {
  
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
  warn <- function(w, warn_txt) {
    warn_msg <- conditionMessage(w)
    warn_con <- file(warn_file, open = "a")
    writeLines(paste(warn_txt, ":", warn_msg), warn_con)
    close(warn_con)
    invokeRestart(findRestart("muffleWarning"))
  }
  
  err <- function(e, err_txt) {
    err_msg <- conditionMessage(e)
    err_con <- file(err_file, open = "a")
    writeLines(paste(err_txt, ":", err_msg), err_con)
    close(err_con)
  }
  
  sp_dirs <- list.dirs(paste0(hv.dir, "/projections/", hv.method))
  
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]
  
  # Load the cavm region
  shapefiles = c(
    cavm = "./resources/region/cavm-noice/cavm-noice.shp"
  )
  
  regions <- import_regions(shapefiles, "./outputs/visualize/region")
  
  cavm <- handle_region(regions$region.cavm)
  
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
  source_all("./src/visualize/components")
  
  visualize_freqpoly(
    sp_cells = freq_stack, 
    region = cavm, 
    region.name = "CAVM", 
    plot.x = "richness",
    plot.color = "country", 
    plot.shade = "floristicProvince"
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
    region.name = "CAVM",
    projection  = laea_crs,
    projection.method = "near"
  )
  
  inc_cover <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    dir = paste0(get_vis_log, "/coverage"), 
    region = NULL,
    n = 9,
    out.order = "coverage",
    fun = get_inclusion_coverage,
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
    region.name = "CAVM",
    extent = cavm_laea_ext,
    projection  = laea_crs,
    projection.method = "near",
    verbose = F
  )
  
  ##########################
  #        Figure 3        #
  ##########################
  
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  # Get the max values of all probability raster files
  prob_max <- parallel_spec_handler(
    spec.dirs = sp_dirs,
    dir = paste0(get_vis_log, "/prob-max"),
    region = NULL,
    hv.project.method  = "probability",
    n = 3,
    out.order = "maxValue",
    fun = get_inclusion_coverage,
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
    region.name = "CAVM"
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
  
  richness_dt <- order_by_apg(richness_dt, by = "order")
  
  print(head(richness_dt, 3))
  
  visualize_richness(
    dt = richness_dt, 
    axis.x = "floristicProvince", 
    axis.y = "relativeRichness", 
    fill = "order",
    group = "order",
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
