main <- function(
    spec.known = NULL,
    spec.unknown = NULL,
    test = NULL,
    column = "scientificName",
    coord.uncertainty = NULL,
    region = NULL,
    download.key = NULL,
    download.doi = NULL,
    hv.iterations = NULL,
    hv.method = "box",
    hv.accuracy = "accurate",
    hv.dims = NULL,
    hv.incl.threshold = 0.5,
    vis.shape = NULL,
    verbose = FALSE
  ) {
  
  if (is.null(spec.unknown) & is.null(test)) {
    vebcat("Error: cannot have both spec.unknown and test as NULL.", color = "fatalError")
    stop("Add function to spec.unknown or use test = 'small' or test = 'big'.")
  }
  
  if (!is.null(test)) {
    if (test == "small") {
      hv_dir <- paste0("./outputs/hypervolume/test-small")
      vis_dir <- paste0("./outputs/visualize/test-small")
    } else if (test == "big") {
      hv_dir <- paste0("./outputs/hypervolume/test-big")
      vis_dir <- paste0("./outputs/visualize/test-big")
    }
  } else {
    hv_dir <- paste0("./outputs/visualize/", gsub("filter_", "", deparse(substitute(spec.unknown))))
    vis_dir <- paste0("./outputs/visualize/", gsub("filter_", "", deparse(substitute(spec.unknown))))
  }
  
  max_cores <- calc_num_cores(
    ram.high = total_cores, 
    verbose =  FALSE
  )
  
  setup_sequence(
    hv.method = "box",
    hv.accuracy = "accurate", 
    hv.incl.t = 0.5,
    hv.dims = c(18, 10, 3, 4),
    cores.max = max_cores$total,
    verbose = TRUE
  )
  
  sp_dir <- filter_sequence(
    spec.known = spec.known, 
    spec.unknown = spec.unknown,
    test = test,
    column = column,
    coord.uncertainty = coord.uncertainty,
    cores.max = max_cores,
    region = region,
    download.key = download.key,
    download.doi = download.doi,
    verbose = verbose
  )
  
  # Get file_names as a list
  catn("Listing files in the directory:", highcat(sp_dir))
  
  sp_list <- list.files(sp_dir, full.names = TRUE)
  
  catn("Found", highcat(length(sp_list)), "species.")
  
  min_disk_space <- get_disk_space("/export", units = "GB") * 0.2
  
  peak_ram <- setup_hv_sequence(
    hv.method = hv.method, 
    hv.accuracy = hv.accuracy, 
    hv.dims = hv.dims, 
    hv.incl.threshold = hv.incl.threshold, 
    verbose = verbose
  )
  
  max_cores <- calc_num_cores(
    ram.high = peak_ram$high, 
    ram.low = peak_ram$low, 
    verbose =  FALSE
  )
  
  cores_max_high <- min(length(sp_list), max_cores$high)
  cores_max_total <- min(length(sp_list), max_cores$total)
  
  vebprint(max_cores$total, text = "Total cores input into the Hypervolume sequence:")
  vebprint(max_cores$high, text = "High load cores input into the Hypervolume sequence:")
  vebprint(max_cores$low, text = "Low load cores input into the Hypervolume sequence:")
  
  # Run the data_acquisition here instead of inside each node.
  hypervolume_sequence(
    spec.list = sp_list,
    iterations = hv.iterations, 
    cores.max.high = cores_max_high,
    cores.max = cores_max_total,
    min.disk.space = min_disk_space,
    verbose = verbose,
    hv.dir = hv_dir,
    hv.method = hv.method, #box approx 13 min, gaussian 1 hours 10 minutes, 
    hv.accuracy = hv.accuracy, 
    hv.dims = hv.dims, 
    hv.incl.threshold = hv.incl.threshold
  )
  
  #as.numeric(gsub("node", "", readLines("outputs/hypervolume/sequence/logs/node-iterations.txt")))
  
  visualize_sequence(
  out.dir = vis_dir,
    shape = vis.shape,
    hv.dir = hv_dir, 
    hv.method = hv.method, 
    projection = "laea",
    verbose = verbose
  )
}
