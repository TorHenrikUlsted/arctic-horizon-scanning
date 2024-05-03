source_all("./src/filter/components")
source_all("./src/filter/custom_filters")

filter_sequence <- function(spec.known = NULL, spec.unknown  = NULL, test = NULL, column = "scientificName", chunk.name = "species", coord.uncertainty =  NULL, cores.max = 1, region = NULL, download.key = NULL, download.doi = NULL, chunk.size = 1e6, chunk.iterations = NULL, verbose = FALSE) {
  
  on.exit(closeAllConnections())
  vebcat("Initiating filtering sequence", color = "seqInit")
  
  ##################################################
  #                     setup                      #
  ##################################################
  
  filter_timer = start_timer("filter_timer")
  
  coord_un_file <- "./outputs/setup/region/coordinateUncertainty-m.txt"
  
  if (!file.exists(coord_un_file)) {
    mdwrite(
      post_seq_nums,
      heading = "1;Filter Sequence"
    )
  }
  
  vebcat("Loading dfs.", veb = verbose)
  
  dts <- select_wfo_column(
    filepath = "./resources/synonym-checked", 
    col.unique = column, 
    col.select = NULL,
    verbose = verbose
  )
  
  dts <- fix_nomatches(
    dfs = dts, 
    nomatch.edited = "./resources/manual-edit/wfo-nomatch-edited.csv", 
    column = column,
    verbose = verbose
  )
  
  if (is.null(coord.uncertainty)) {
    coord.uncertainty <- as.numeric(readLines(coord_un_file))
  }
    
    ####################
    # Filter lists  
    ####################
  
  if (!is.null(test)) {
    spec.known <- get("filter_test_known")
    
    if (test == "small") {
      spec.unknown <- get("filter_test_small")
      download.key = "0180552-240321170329656"
      download.doi = "https://doi.org/10.15468/dl.xzxdpx"
    } else if (test == "big") {
      spec.unknown <- get("filter_test_big")
      download.key = "0180885-240321170329656"
      download.doi = "https://doi.org/10.15468/dl.sgf54g"
    } else {
      vebcat("Test has to be either 'small' or 'big'.", color = "fatalError")
      stop("Change the test parameter.")
    }
  } 
  print(is.null(spec.known))
  if (!is.null(spec.known)) {
    known <- spec.known(
      dts = dts,
      column = column,
      verbose = verbose
    )
  } else {
    known <- NULL
  }
  
  vebprint(known, text = "known output:", veb = verbose)
  
  if (is.null(spec.unknown)) {
    vebcat("Error, spec.unknown is NULL.", color = "fatalError")
    return(NULL)
  } else {
    unknown <- spec.unknown(
      known.filtered = known,
      dts = dts,
      column = column,
      verbose = verbose
    )
  }
  
    ###############################
    # Unknown Occurrence download
    ###############################
  
  occ_name <- paste0(unknown$dir, "/", basename(unknown$dir), "-occ")
  
  unknown_occ <- get_occurrence(
    spec = unknown$spec,
    file.out = occ_name,
    region = region,
    coord.uncertainty = coord.uncertainty,
    download.key = download.key,
    download.doi = download.doi,
    verbose = verbose
  )
  
  if (!is.null(unknown_occ$occ)) {
    sp_occ <- unknown_occ$occ
  } else {
    sp_occ <- paste0(occ_name, ".csv")
  } 
  
  #####################
  # Chunking protocol 
  #####################
  
  chunk_dir <- paste0(unknown$dir, "/chunk")
  
  chunk_protocol(
    spec.occ = sp_occ,
    spec.keys = unknown_occ$keys,
    chunk.name = chunk.name,
    chunk.col = "cleanName",
    chunk.dir = chunk_dir,
    chunk.size = chunk.size,
    cores.max = 1,
    iterations = chunk.iterations,
    verbose = verbose
  )
  
  end_timer(filter_timer)
  
  vebcat("Filtering sequence completed successfully.", color = "seqSuccess")
  
  return(paste0(chunk_dir, "/", chunk.name))
}
