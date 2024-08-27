check_syn_wfo <- function(checklist, column, out.dir, cores.max = 1, verbose = FALSE, counter = 1) {
  if (!"data.table" %in% class(checklist) && !"data.frame" %in% class(checklist)) {
    stop("The input data is not in the 'data.table' or 'data.frame' format.", print(class(checklist)))
  }
  checklist <- as.data.table(checklist)
  out_file <- paste0(out.dir, "/wfo-match.csv")
  
  if (file.exists(out_file)) {
    catn("Reading existing WFO.match file from:", colcat(out_file, color = "output"))
    return(fread(out_file))
  }
  
  vebcat("Initiating WFO synonym check", color = "funInit")
  wfo_timer <- start_timer("wfo_match")
  
  catn("Running the WFO synonym check for", highcat(nrow(checklist)), "species with column", highcat(column))
  vebprint(head(checklist, 3), text = "Table sample:")
  
  if (nrow(checklist) < 10) {
    wfo_result <- WFO.match(
      spec.data = checklist, 
      spec.name = column, 
      WFO.file = WFO_file, 
      verbose = verbose, 
      counter = counter
    )
  } else {
    
    custom_evals <- list(
      packages = c(
        "WorldFlora", 
        "data.table"
      ),
      source = c(
        "./src/utils/components/condition_handlers.R",
        "./src/utils/components/time_tracker.R",
        "./src/utils/components/file_managers.R",
        "./src/utils/components/loader.R"
      )
    )
    
    cores_max <- calc_num_cores(
      ram.high = 2,
      cores.total = config$memory$total_cores,
      verbose = verbose
    )
    
    wfo_result <- wfo_parallel(
      checklist = checklist,
      column = column,
      out.dir = out.dir,
      cores.max = min(nrow(checklist), cores_max$total),
      evals = custom_evals,
      counter = counter,
      verbose = verbose
    )
  }
  
  
  if (!is.data.table(wfo_result)) wfo_result <- as.data.table(wfo_result)
  fwrite(wfo_result, paste0(out.dir, "/wfo-match.csv"), bom = T)
  
  wfo_result <- wfo_mismatch_check(
    wfo.result = wfo_result, 
    col.origin = column,
    out.file = paste0(out.dir, "/wfo-match-mismatches.csv"),
    unchecked = FALSE,
    verbose = verbose
  ) # returns list of clean and mismatched info
  
  end_timer(wfo_timer)
  
  vebcat("WFO.match synonym check completed", color = "funSuccess")
  
  return(wfo_result)
}
