main <- function(
    spec.known = NULL,
    spec.known.key = NULL,
    spec.known.doi = NULL,
    spec.unknown = NULL,
    spec.unknown.key = NULL,
    spec.unknown.doi = NULL,
    gbif.occ.region = NULL,
    coord.uncertainty = NULL,
    hv.iterations = NULL,
    hv.method = "box",
    hv.accuracy = "accurate",
    hv.dims = NULL,
    hv.incl.threshold = 0.5,
    vis.shape = NULL,
    vis.projection = "laea",
    vis.title = TRUE,
    vis.region.name = "Region", 
    vis.subregion.name = "Sub Region", 
    vis.composition.taxon = "order", 
    vis.save.device = "jpeg",
    vis.save.unit = "px",
    plot.show = FALSE,
    validation = TRUE,
    total.cores = 1,
    verbose = FALSE,
    force.seq = NULL
  ) {
  
  run <- 1
  
  repeat {
    if (is.null(force.seq) && run == 1 || force.seq != "validation" && run == 1) {
      validation_run <- FALSE
    } else if (!is.null(force.seq) && force.seq == "validation" || run == 2) {
      validation_run <- TRUE
      vebcat("Initiating Validation Protocol", color = "funInit")
    } else {
      vebcat("Error when trying to set validation parameter in main, stopping...", color = "fatalError")
      vebprint(run, text = "run:")
      vebprint(force.seq, text = "force.seq:")
      stop("Run more than two, stopping indefinte loop.")
    }
    
    if (!validation_run) {
      hv_dir <- paste0("./outputs/hypervolume/", spec.unknown)
      vis_dir <- paste0("./outputs/visualize/", spec.unknown)
    } else {
      hv_dir <- paste0("./outputs/hypervolume/", paste0(spec.unknown, "_validation"))
      vis_dir <- paste0("./outputs/visualize/", paste0(spec.unknown, "_validation"))
    }
      
    vis.shape = paste0("./outputs/setup/region/", vis.shape, "/", vis.shape, ".shp")
    
    max_cores <- calc_num_cores(
      ram.high = total.cores,
      cores.total = total.cores,
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
      spec.known = list(
        name = spec.known,
        download.key = spec.known.key,
        download.doi = spec.known.doi
      ),
      spec.unknown = list(
        name = spec.unknown,
        download.key = spec.unknown.key,
        download.doi = spec.unknown.doi
      ),
      validation = validation_run,
      coord.uncertainty = coord.uncertainty,
      cores.max = max_cores,
      region = gbif.occ.region,
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
      cores.total = total.cores,
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
      vis.projection = vis.projection,
      vis.title = vis.title,
      vis.region.name = vis.region.name,
      vis.subregion.name = vis.subregion.name,
      vis.composition.taxon = vis.composition.taxon,
      vis.save.device = vis.save.device,
      vis.save.unit = vis.save.unit,
      plot.show = plot.show,
      verbose = verbose
    )
    
    rm(sp_dir, sp_list, min_disk_space, peak_ram, max_cores, cores_max_high, cores_max_total)
    invisible(gc())
    
    if (run == 2 || !is.null(force.seq) && force.seq == "validation") {
      catn("Validation analysis has finished, closing main loop.\n")
      break
    }
    
    if (validation) {
      run <- 2
      catn("Changing mode to validation")
    } else {
      catn("Validation analysis is set to FALSE, closing main loop.\n")
      break
    }
  }
}
