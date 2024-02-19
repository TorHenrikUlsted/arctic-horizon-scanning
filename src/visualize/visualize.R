source_all("./src/visualize/components")

visualize <- function(spec.list, out.dir, hv.dir, hv.method, x.threshold, verbose) {
  
  # Test params
  spec.list = list.files("./outputs/filter/arctic/chunk/species", full.names = TRUE)
  out.dir = "./outputs/visualize" 
  hv.dir = "./outputs/hypervolume/sequence"
  hv.method = "box"
  hv.projection = laea_crs
  x.threshold = 0.2
  show.plot = FALSE
  verbose = T
  
  # Set up directories
  log_dir <- paste0(out.dir, "/logs")
  create_dir_if(log_dir)
  
  # Set up error handling
  warn_file <- paste0(log_dir, "/", hv.method, "-warning.txt")
  err_file <- paste0(log_dir, "/", hv.method, "-error.txt")
  
  create_file_if(warn_file)
  create_file_if(err_file)
  
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
  
  ##########################
  #       Check data       #
  ##########################
  source_all("./src/visualize/components")
  # Check the output if it is in the expected format and length
  checked_data <- check_hv_output(spec.list, hv.dir, hv.method, hv.projection, hv.inc.t = 0.5, verbose = F)
  
  # Clean the output data - 9 hours
  cleaned_data <- clean_hv_output(checked_data, out.dir, hv.dir, hv.projection, verbose = F)
  
  # Make excluded species list
  visualize_data <- get_visualize_data(cleaned_data, hv.dir, hv.method, verbose = F, warn = warn, err = err)

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
  
  ### THE APP() GOES INSANELY MUCH SLOWER WITH LAEA THAN LONGLAT
  
  nlyr(visualize_data$inc_stack)
  
  test_rast <- subset(visualize_data$inc_stack, 1:10)
  
  test_rast <- terra::project(test_rast, hv.projection)
  
  test_prob <- subset(visualize_data$prob_stack, 1:10)
  
  test_prob <- terra::project(test_prob, hv.projection)
  
  test_cavm <- project(cavm, hv.projection)
  
  # Calc app ram usage
  app_ram_peak <- setup_app_process(cavm)
  
  cmu <- get_mem_usage(type = "total", format = "gb") * 0.7
  terraOptions(memfrac = 0.8)
  
  cores_inc_max <- floor(cmu / app_ram_peak$peak_mem_inc)
  cores_prob_max <- floor(cmu / app_ram_peak$peak_mem_prob)
  
  inc_timer <- start_timer("inc_cell_timer")
  inc_cells <- get_sp_cell(test_rast, test_cavm, method = "inclusion", cores = 1, out.filename = paste0(out.dir, "/rast/test-inc-laea.tif"), verbose = T)
  end_timer(inc_timer)
  
  prob_timer <- start_timer("Prob_timer")
  prob_cells <- get_sp_cell(test_prob, test_cavm, method = "probability", cores = 1, out.filename = paste0(out.dir, "/rast/test-prob.tif"))
  end_timer(prob_timer)
  
  # Create figure 1
  # Reproject the stacks
  laea_inc_timer <- start_timer("Laea inc")
  inc_cells$sp_count <- project(inc_cells$sp_count, laea_crs)
  end_timer(laea_inc_timer)
  
  cavm_timer <- start_timer("reproject cavm")
  cavm <- project(cavm, laea_crs)
  end_timer(cavm_timer)
  
  laea_prob_timer <- start_timer("Laea prob")
  prob_cells$sp_count <- project(prob_cells$sp_count, laea_crs)
  end_timer(laea_prob_timer)
  
  # Get cell sums
  test_rast <- ceiling(test_rast)
  
  # Get sub_region sums
  t_ex <- terra::extract(test_rast, cavm, cells = TRUE)
  t_ex <- as.data.table(t_ex, na.rm = TRUE)
  
  # Add cavm regions
  cavm_dt <- as.data.table(cavm)
  cavm_dt[, ID := .I]
  
  sub_regions_dt <- merge(t_ex, cavm_dt[, .(ID, FLOREG, country, floristicProvince)], by = "ID")
  
  # Calculate the sum of species in each cell and in each cell in each floreg
  sp_cols <- 2:(1 + nlyr(test_rast))
  
  sub_regions_dt <- sub_regions_dt[, cellSum := rowSums(.SD, na.rm = TRUE), .SDcols = sp_cols]
  
  #sub_regions_dt <- sub_regions_dt[, floregSum := sum(cellSum, na.rm = TRUE), by = .(FLOREG)]
  
  #t_ex_long <- merge(t_ex_long, stat_inc_data, by = "species", all.x = TRUE, allow.cartesian = TRUE)
  
  visualize_histogram(sp_cells = sub_regions_dt, region = cavm, region.sub.color = "country", region.sub.shades = "floristicProvince", region.name = "CAVM")
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  visualize_hotspots(sp_cells = final_dt, prob.value = "proportion", region = cavm, region.sub = "country", region.name = "CAVM", plot.show = T) 
  
  cavm_sf <- st_as_sf(cavm)
  final_sf <- merge(cavm_sf, final_dt, by.x = "FLOREG", by.y = "FLOREG")
  
  ggplot() +
    geom_sf(data = final_sf, aes(fill = cell_sum)) +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, guide = guide_legend(reverse = TRUE)) +
    labs(x = "Longitude", y = "Latitude", title = paste0("Potential species distribution in the ", region.name), fill = "Proportion", color = "Country") +
    coord_sf(xlim = c(region_ext$xmin, region_ext$xmax), ylim = c(region_ext$ymin, region_ext$ymax)) +
    theme_minimal() +
    labs(fill = "Cell Sum",
         title = "Hotspots on CAVM Shape",
         caption = "Red areas indicate higher cell sums")
  

  # Figure 4: Matrix with floristic regions 
  visualize_matrix(sp_cells = inc_cells, region.sub = "country")
  

  
  source_all("./src/visualize/components")
  # Figure 5: Sankey with floristic regions 
  visualize_sankey(dt.src = visualize_data$included_sp, dt.target = sub_regions_dt)
  
}
