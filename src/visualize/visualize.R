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
  
  create_file_if(warn_file)
  create_file_if(err_file)
  
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
  
  cavm <- regions$cavm
  
  cavm_desc <- fread("./resources/region/cavm2003-desc.csv")
  
  index <- match(cavm$FLOREG, cavm_desc$FLOREG)
  
  cavm$country <- cavm_desc$country[index]
  cavm$floristicProvince <- cavm_desc$floristicProvince[index]
  
  # Remove the ice sheet
  cavm <- cavm[cavm$FLOREG != 0, ]
  cavm <- na.omit(cavm)
  
  cavm_laea <- terra::project(cavm, laea_crs)
  
  cavm_laea_ext <- ext(cavm_laea)
  
  ##########################
  #       Check data       #
  ##########################
  
  # Check the output if it is in the expected format and length
  checked_data <- check_hv_output(spec.list, hv.dir, hv.method, hv.projection, hv.inc.t = 0.5, verbose = F)
  
  # Clean the output data - 9 hours
  cleaned_data <- clean_hv_output(checked_data, out.dir, hv.dir, hv.projection, verbose = T)
  
  # Get stats csv file
  sp_stats <- get_stats_data(cleaned_data, hv.dir, hv.method, verbose = F, warn = warn, err = err)
  
  # Get the max values of all probability raster files
  prob_max_values <- get_projection_max(hv.dir, hv.method, hv.project.method = "probability", n = 9)
  
  # Read the rasters with the highest max values
  prob_stack <-  stack_projections(filenames = prob_max_values, projection = laea_crs, projection.method = "bilinear")
  # Make excluded species list
  inc_stack <- get_inclusion_data(region = cavm, hv.dir, hv.method, hv.project.method = "inclusion-0.5", verbose = F, warn = warn, err = err)
  
  # Get number of cells per floristic region
  freq_stack <- na.omit(inc_stack)
  
  freq_stack <- freq_stack[, ncellsRegion := .N, by = floristicProvince]
  
  freq_stack <- freq_stack[, propCells := ncellsRegion / nrow(freq_stack), by = floristicProvince]
  
  visualize_freqpoly(
    sp_cells = freq_stack, 
    region = cavm, 
    region.name = "CAVM", 
    plot.x = "richness",
    plot.color = "country", 
    plot.shade = "floristicProvince"
    )
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  
  # Convert to raster by using a template
  hotspot_raster <- convert_template_raster(
    input.values = inc_stack$richness,
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
  
  inc_cover_t <- parallel_spec_handler(
    spec.dirs = sp_dirs, 
    sub.dir = paste0(get_vis_log, "/coverage"), 
    n = 3,
    fun = get_inclusion_coverage,
    verbose = FALSE
  )
  
  inc_cover <- stack_projections(
    filenames = inc_cover, 
    projection = laea_crs, 
    projection.method = "near", 
    binary = TRUE
  )
  
  visualize_highest_spread(
    rast = inc_cover,
    region = world_map,
    region.name = "CAVM",
    extent = cavm_laea_ext,
    projection  = laea_crs,
    projection.method = "near",
    verbose = F
  )
  
  source_all("./src/visualize/components")
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  prob_max <- parallel_spec_handler(
    spec.dirs = sp_dirs,
    hv.project.method  = "probability",
    sub.dir = paste0(get_vis_log, "/prob-max"), 
    n = 9,
    fun = get_inclusion_coverage,
    verbose = FALSE
  )
  
  prob_stack <- stack_projections(
    filenames = prob_max, 
    projection = laea_crs, 
    projection.method = "near", 
    binary = TRUE
  )
  
  visualize_suitability(
    rast = prob_stack, 
    region = cavm_laea, 
    region.name = "CAVM"
  )

  # Figure 4: stacked barplot wtih taxa per floristic region
  sp_regions_dt <- get_sp_names_region()
  
  richness_dt <- calculate_richness(sub_regions_dt, sp_cols, taxon = "family")
  
  print(head(richness_dt, 3))
  
  visualize_richness(
    dt = richness_dt, 
    axis.x = "floristicProvince", 
    axis.y = "relativeRichness", 
    fill = "family",
    group = "family"
  )
  
  # Figure 5: Sankey with floristic regions 
  visualize_sankey(
    dt.src = visualize_data$included_sp, 
    dt.target = sub_regions_long
    )
  
}
