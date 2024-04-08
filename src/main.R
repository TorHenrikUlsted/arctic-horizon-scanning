main <- function(
    spec.known = NULL,
    speck.unknown,
    test = NULL,
    column = "scientificName",
    coord.uncertainty = 4600,
    region = NULL,
    download.key = NULL,
    download.doi = NULL,
    hv.iterations = NULL,
    hv.method = "box",
    hv.accuracy = "accurate",
    hv.dims = NULL,
    hv.incl.threshold = 0.5,
    verbose = FALSE
  ) {
  
  
  max_cores <- calc_num_cores(
    ram.high = total_cores, 
    verbose =  FALSE
  )
  
  setup_sequence(
    hv.method = "box",
    hv.accuracy = "accurate", 
    hv.incl.t = 0.5,
    hv.dims = c(18, 10, 3, 4),
    cores.max = max_cores,
    verbose = FALSE
  )
  
  sp_dir <- filter_sequence(
    # The function used to get species known in the region
    spec.known = filter_arctictest, 
    # this function uses the spec.known to remove from spec.absent
    spec.unknown = filter_glonaftest,
    test = "small",
    column = "scientificName",
    coord.uncertainty = as.numeric(readLines("./outputs/hypervolume/data_acquisition/logs/coordinateUncertainty-m.txt")),
    cores.max = max_cores,
    region = NULL,
    download.key = NULL,
    download.doi = NULL,
    verbose = FALSE
  )
  
  # Get file_names as a list
  catn("Listing files in the directory:", highcat(sp_dir))
  
  sp_list <- list.files(sp_dir, full.names = TRUE)
  
  catn("Found", highcat(length(sp_list)), "species.")
  
  min_disk_space <- get_disk_space("/export", units = "GB") * 0.2
  
  # Check memory peak of one node by conducting a test run as well as setting up the entire hypervolume sequence
  peak_ram <- setup_hv_sequence(min_disk_space)
  
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
    iterations = NULL, 
    cores.max.high = cores_max_high,
    cores.max = cores_max_total,
    min.disk.space = min_disk_space,
    verbose = TRUE,
    hv.method = "box", #box approx 13 min, gaussian 1 hours 10 minutes, 
    hv.accuracy = "accurate", 
    hv.dims = c(18, 10, 3, 4), 
    hv.incl.threshold = 0.5
  )
  
  #as.numeric(gsub("node", "", readLines("outputs/hypervolume/sequence/logs/node-iterations.txt")))
  
  visualize_sequence(
    out.dir = "./outputs/visualize", 
    hv.dir = "./outputs/hypervolume/sequence", 
    hv.method = "box", 
    verbose = FALSE
  )
}
