source_all("./src/visualize/components")

visualize <- function(out.dir, hv.dir, hv.method, x.threshold, verbose) {
  
  # Test params
  out.dir = "./outputs/visualize" 
  hv.dir = "./outputs/hypervolume/sequence"
  hv.method = "box"
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

  visualize_data <- get_visualize_data(hv.dir, hv.method, verbose = T, warn = warn, err = err)
  
  # Load regions
  
  shapefiles = c(
    cavm = "./resources/region/cavm-noice/cavm-noice.shp"
  )
  
  regions <- import_regions(shapefiles, "./outputs/visualize/region")
  
  # cavm_floreg <- terra::split(regions$cavm, regions$cavm$VEGPHYS)

   # for (i in seq_along(cavm_floreg)) {
   #   names(cavm_floreg)[i] <- paste("floreg", unique(cavm_floreg[[i]]$FLOREG), sep = "_")
   #   if (verbose) cat("Renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(names(cavm_floreg)[i]), "\n")
   # }
  
  regions$cavm
  plot(test_rast[[1]])
  plot(regions$cavm)
  plot(cavm_floreg$floreg_0)
  
  test_rast[[1]]
  cavm_floreg[[1]]
  
  # Check raster files for different things and crs
  
  nlyr(visualize_data$inc_stack)
  
  test_rast <- subset(visualize_data$inc_stack, 1:10)
  
  
  terra::cellSize(test_rast, unit="km", lyrs=TRUE, mask=TRUE)
  
  terra::res(test_rast)
  
  prod(terra::res(test_rast))
  
  # Get species per cell
  
  sp_cell <- get_sp_cell()
  
  source_all("./src/visualize/components")
  
  # Maybe use get_sp_cell
  
  # Create figure 1
  
  # Missing :: Add floristic region names, and combine the ones in the different regions to countries
  make_histogram(rast = test_rast, region = regions$cavm, region.sub = "VEGPHYS", region.name = "CAVM")
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  
  
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  
  # Figure 4: Matrix with floristic regions 
  
  # Figure 5
  
}