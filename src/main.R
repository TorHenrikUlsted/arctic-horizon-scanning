main <- function(
    spec.known = NULL,
    spec.unknown = NULL,
    coord.uncertainty = NULL,
    gbif.occ.region = NULL,
    download.key = NULL,
    download.doi = NULL,
    hv.iterations = NULL,
    hv.method = "box",
    hv.accuracy = "accurate",
    hv.dims = NULL,
    hv.incl.threshold = 0.5,
    vis.shape = NULL,
    vis.title = TRUE,
    vis.region.name = "Region", 
    vis.subregion.name = "Sub Region", 
    vis.composition.taxon = "order", 
    vis.save.device = "jpeg",
    vis.save.unit = "px",
    plot.show = FALSE,
    verbose = FALSE,
    force.seq = NULL
  ) {
  
  hv_dir <- paste0("./outputs/hypervolume/", spec.unknown)
  vis_dir <- paste0("./outputs/visualize/", spec.unknown)
  vis.shape = paste0("./outputs/setup/region/", vis.shape, "/", vis.shape, ".shp")
  
  max_cores <- calc_num_cores(
    ram.high = config$memory$total_cores,
    cores.total = config$memory$total_cores,
    verbose = verbose
  )
  
  setup_sequence(
    hv.method = hv.method,
    hv.accuracy = hv.accuracy, 
    hv.incl.threshold = hv.incl.threshold,
    hv.dims = hv.dims,
    cores.max = max_cores$total,
    force.seq = force.seq,
    verbose = verbose
  )
  
  if (is.null(spec.unknown)) {
    vebcat("Error: cannot have spec.unknown as NULL", color = "fatalError")
    stop("Add name to spec.unknown or use 'test_small' OR 'test_big'.")
  }
  
  sp_dir <- filter_sequence(
    spec.known = spec.known, 
    spec.unknown = spec.unknown,
    coord.uncertainty = coord.uncertainty,
    cores.max = max_cores,
    region = gbif.occ.region,
    download.key = download.key,
    download.doi = download.doi,
    force.seq = force.seq,
    verbose = verbose
  )
  
  # Get file_names as a list
  catn("Listing files in the directory:", highcat(sp_dir))
  
  sp_list <- list.files(sp_dir, full.names = TRUE)
  
  catn("Found", highcat(length(sp_list)), "species.")
  
  min_disk_space <- get_disk_space("/export", units = "GB") * 0.2
  
  invisible(gc())
  
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
    cores.total = config$memory$total_cores,
    verbose =  verbose
  )
  
  cores_max_high <- min(length(sp_list), max_cores$high)
  cores_max_total <- min(length(sp_list), max_cores$total)
  
  vebcat("Total cores input into the Hypervolume sequence:", highcat(max_cores$total))
  vebcat("High load cores input into the Hypervolume sequence:", highcat(max_cores$high))
  vebcat("Low load cores input into the Hypervolume sequence:", highcat(max_cores$low))
  
  # Run the data_acquisition here instead of inside each node.
  hypervolume_sequence(
    spec.list = sp_list,
    iterations = hv.iterations, 
    cores.max.high = cores_max_high,
    cores.max = cores_max_total,
    min.disk.space = min_disk_space,
    verbose = verbose,
    hv.dir = hv_dir,
    hv.method = hv.method, #box approx 20 min, gaussian 1 hour 15 minutes, 
    hv.accuracy = hv.accuracy, 
    hv.dims = hv.dims, 
    hv.incl.threshold = hv.incl.threshold
  )
  
  visualize_sequence(
    out.dir = vis_dir,
    res.unknown = spec.unknown,
    res.known = spec.known,
    shape = vis.shape,
    hv.dir = hv_dir, 
    hv.method = hv.method,
    vis.projection = config$simulation$projection,
    vis.title = vis.title,
    vis.region.name = vis.region.name,
    vis.subregion.name = vis.subregion.name,
    vis.composition.taxon = vis.composition.taxon,
    vis.save.device = vis.save.device,
    vis.save.unit = vis.save.unit,
    plot.show = plot.show,
    verbose = verbose
  )
  
  
}
