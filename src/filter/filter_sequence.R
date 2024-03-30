source_all("./src/filter/components")
source_all("./src/filter/custom_filters")

filter_sequence <- function(spec.known, spec.unknown, column = "scientificName", chunk.name = "species", test = NULL, cores.max = 1, verbose = FALSE) {
  on.exit(closeAllConnections())
  vebcat("Initiating filtering sequence", color = "seqInit")
  
  
  ##################################################
  #                     setup                      #
  ##################################################
  
  filter_timer = start_timer("filter_timer")
  on.exit(end_timer(filter_timer))
  
  vebcat("Loading dfs.", veb = verbose)
  
  dfs <- select_wfo_column(
    filepath = "./resources/synonym-checked", 
    col.unique = column, 
    col.select = NULL,
    verbose = T
  )
  
  dfs <- fix_nomatches(
    dfs = dfs, 
    nomtach.edited = "./resources/manual-edit/wfo-nomatch-edited.csv", 
    column = column
  )
    
    ####################
    # Filter lists  
    ####################
    
  if (!is.null(spec.known)) {
    known <- spec.known(
      dfs = dfs, 
      verbose = verbose
    )
  } else {
    known_species <- NULL
  }
  
  unknown <- spec.unknown(
    spec.known = known,
    dfs = dfs, # This is if you have multiple dfs and need to remove to avoid duplicates
    verbose = verbose
  )
    
    ###############################
    # Unknown Occurrence download
    ###############################
  
  occ_name <- paste0(unknown$occ.dir, "/", basename(unknown$occ.dir), "-occ")

  unknown_occ <- get_occurence(
    spec = unknown_species$spec,
    file.out = occ_name,
    region = NULL,
    download.key = NULL,
    download.doi = NULL,
    verbose = verbose
  )
  
  if (!is.null(unknown_occ)) {
    sp_occ <- unknown_occ
  } else {
    sp_occ <- paste0(occ_name, ".csv")
  } 
    
  
  #####################
  # Chunking protocol 
  #####################
  
  chunk_dir <- paste0(unknown$dir, "/chunk")
  
  chunk_protocol(
    spec.occ = sp_occ,
    chunk.name = chunk.name,
    chunk.col = column,
    chunk.dir = chunk_dir,
    chunk.size = 1e6,
    cores.max = 1,
    iterations = NULL,
    verbose
  )
  
  vebcat("Filtering sequence completed successfully.", color = "seqSuccess")
  
  return(paste0(chunk_dir, "/", chunk.name))
}
