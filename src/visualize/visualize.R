source_all("./src/visualize/components")

visualize <- function(spec.list, out.dir, hv.dir, hv.method, x.threshold, verbose) {
  
  # Test params
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
  
  # Make excluded species list

  visualize_data <- get_visualize_data(hv.dir, hv.method, verbose = F, warn = warn, err = err)
  
  ## ADD A CHECK BEFORE RETURNING DATA IN THE HV ANALYSIS
  
  ##########################
  #       Check data       #
  ##########################
  source_all("./src/visualize/components")
  
  checked_data <- check_hv_output(spec.list, hv.dir, hv.method, hv.projection, hv.inc.t = 0.5, hv.clean = F, verbose = T)
  
  cleaned_data <- clean_hv_output(visualize_data, projection = laea_crs, verbose = T)
  
  # Load regions
  vals <- c(0, 1, 2, 3, 4, 5, 6, 7, 81, 82, 83, 84, 91, 92, 93, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)
  
  shapefiles = c(
    cavm = "./resources/region/cavm-noice/cavm-noice.shp"
  )
  
  regions <- import_regions(shapefiles, "./outputs/visualize/region")
  
  cavm_floreg <- terra::split(regions$cavm, regions$cavm$FLOREG)
  cavm_floreg <- svc(cavm_floreg)
  
  cavm <- terra::project(regions$cavm, laea_crs)
  
  floreg <- terra::split(cavm, cavm$FLOREG)
  
  for (i in seq_along(cavm_floreg)) {
    names(cavm_floreg)[i] <- paste("floreg", unique(cavm_floreg[[i]]$FLOREG), sep = "_")
    if (verbose) cat("Renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(names(cavm_floreg)[i]), "\n")
    #plot(cavm_floreg[[i]])
  }
  
  ### NEED TO HANDLE GREENLAND 0 VALUE - check if raster has things we can use
  
  test_rast[[1]]
  
  nlyr(visualize_data$inc_stack)
  
  test_rast <- subset(visualize_data$inc_stack, 1:10)
  test_prob <- subset(visualize_data$prob_stack, 1:10)
  
  # Calc app ram usage
  app_ram_peak <- setup_app_process(cavm)
  
  cmu <- get_mem_usage(type = "total", format = "gb") * 0.8
  
  cores_inc_max <- floor(cmu / app_ram_peak$peak_mem_inc)
  cores_prob_max <- floor(cmu / app_ram_peak$peak_mem_prob)
  
  source_all("./src/visualize/components")
  inc_cells <- get_sp_cell(test_rast, cavm, cores = cores_inc_max, method = "inclusion")
  
  prob_cells <- get_sp_cell(test_prob, cavm, cores = cores_prob_max, method = "probability")
  
  # Create figure 1
  
  # Missing :: Add floristic region names, and combine the ones in the different regions to countries
  
  test_frq <- freq(test_rast, digits=0, value=NULL, bylayer=TRUE, usenames=FALSE, zones=regions$cavm, wide=FALSE)
  
  visualize_histogram(rast = test_rast, region = regions$cavm, region.sub = "FLOREG", region.name = "CAVM")
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  visualize_hotsposts() 
  
  plot(terra::ext(regions$cavm), border = NA, col = NA)
  for (i in 1:length(cavm_floreg)) {
    plot(cavm_floreg[[i]], border = region_palette[i], col = NA, lwd = 0.1, add = TRUE)
  }
  plot(inc_cells$sp_count[['prop']], main="Species Hotspots", col = c(NA, terrain.colors(99)), add = T)
  
  inc_cells$sp_count[inc_cells$sp_count == 0] <- NA
  min_inc <- min(terra::values(inc_cells$sp_count$prop), na.rm = TRUE)
  max_inc <- max(terra::values(inc_cells$sp_count$prop), na.rm = TRUE)
  region_ext <- terra::ext(cavm)
  regions_factor <- 
  
  inc_cells$sp_count <- project(inc_cells$sp_count, laea_crs)
  
  
  fig2 <- ggplot() +
    geom_spatvector(data = cavm, aes(color = as.factor(FLOREG)), fill = NA) +
    scale_color_grey( guide = guide_legend(reverse = TRUE)) +
    
    geom_spatraster(data = inc_cells$sp_count, aes(fill = prop)) +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, breaks = c(seq(min_inc, max_inc, by = 0.2), max_inc), limits = c(min_inc, max_inc), guide = guide_legend(reverse = TRUE)) +
    
    labs(x = "Longitude", y = "Latitude", title = paste0("Potential species distribution in the CAVM"), fill = "Proportion") +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  ggsave("./outputs/visualize/plots/figure-2.png", plot = fig2)
  
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  prob_min <- min(terra::values(test_prob), na.rm = TRUE)
  prob_max <- max(terra::values(test_prob), na.rm = TRUE)
  
  plot(terra::ext(cavm), border = NA, col = NA)
  for (i in 1:length(floreg)) {
    plot(floreg[[i]], border = region_palette[i], col = NA, lwd = 0.1, add = TRUE)
  }
  plot(prob_cells$sp_count[['prob_mean']], main="Mean probabillity of species in each cell", col = c(NA, heat.colors(99)), add = T)
  
  ggplot() +
    geom_spatvector(data = cavm) 
  
  prob_min <- min(terra::values(prob_cells$sp_count$prob_mean), na.rm = TRUE)
  prob_max <- max(terra::values(prob_cells$sp_count$prob_mean), na.rm = TRUE)

  # FIX EXTENT
  fig3 <- ggplot() +
    geom_spatvector(data = cavm) +
    geom_spatraster(data = prob_cells$sp_count, aes(fill = prob_mean)) + #maxcell = ncell(prob_cells$sp_count)
    scale_fill_gradientn(colors = c(NA, heat.colors(100)), na.value = NA, breaks = c(seq(prob_min, prob_max, by = 0.2), prob_max),  limits = c(prob_min, prob_max) ) +
    labs(x = "", y = "Latitude", title = paste0("Predicted species distribution in the CAVM"), fill = "Probability") +
    scale_x_continuous(limits = c(0, NA)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  ggsave("./outputs/visualize/plots/figure-3.png", plot = fig3)
  
  
  # Figure 4: Matrix with floristic regions 
  
  density(prob_cells$sp_count, maxcells=100000, plot=TRUE)
  
  
  
  fig4_mat <- matrix()
  
  
  image(fig4_mat, "Matrix Fgiure", col = prob_palette)
  
  
  # Figure 5 -Sankey
  
}
