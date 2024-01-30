source_all("./src/visualize/components")

visualize <- function(out.dir, hv.dir, hv.method, x.threshold, verbose) {
  
  # Test params
  out.dir = "./outputs/visualize" 
  hv.dir = "./outputs/hypervolume/sequence"
  hv.method = "box"
  x.threshold = 0.2 
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
  
  plot(visualize_data$prob_stack[[2]])
  
  # Figure 1A: Histogram with proportion of overlap statistics
  
  included_sp$thresholdLabels <- ifelse(included_sp$overlapRegion > 0.2, included_sp$species, "")
  
  lv <- setNames(included_sp$thresholdLabels, included_sp$species)
  
  ggplot(included_sp, aes(x = species, y = overlapRegion)) +
    geom_col() +
    scale_x_discrete(labels = lv, guide = guide_axis(n.dodge = 2)) +
    geom_label_repel(label = lv) +
    labs(x = "Species", y = "Region overlap") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Figure 1B: Histogram with proportion of overlap statistics and floristic regions
  
  # Figure 2: Stack inclusion tif files and calculate species in each cell to get potential hotspots
  
  # Figure 3: Stack probability tif files and use the highest numbers to get a color gradient
  
  # Figure 4: Matrix with floristic regions 
  
  # Figure 5
  
}