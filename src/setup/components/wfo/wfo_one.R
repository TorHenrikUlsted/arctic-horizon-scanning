check_syn_wfo_one <- function(wfo.match.dt, column = "scientificName", out.dir = "./", verbose = FALSE) {
  vebcat("Initiating the WFO.one synonym check", color = "funInit")
  
  if (!is.data.table(wfo.match.dt)) {
    vebcat("Input data not a data.table", color = "fatalError")
    vebprint(wfo.match.dt, text = "found:")
    stop("Check WFO.one input data")
  }
  
  wfo_one_checklist <- WFO.one(
    WFO.result = wfo.match.dt, 
    priority = "Accepted",
    spec.name = "scientificName",
    verbose = FALSE, 
    counter = 10000
  )
  
  fwrite(wfo_one_checklist, paste0(out.dir, "/wfo-one.csv"), bom = T)
  out_file <- paste0(out.dir, "/wfo-one-clean.csv")
  
  wfo_one_mis <- wfo_mismatch_check(
    wfo.result = wfo_one_checklist, 
    col.origin = column,
    out.file = paste0(out.dir, "/wfo-one-mismatch.csv"),
    unchecked = FALSE,
    verbose = verbose
  ) # returns list of clean and mismatched info
  
  wfo_one_nomatch <- wfo_nomatch_check(
    wfo.result = wfo_one_mis$clean, 
    out.file = paste0(out.dir, "/wfo-one-nomatch.csv"),
    verbose = verbose
  ) # removes nomatches
  
  wfo_one_mis$clean = NULL # remove the clean to save memory
  
  wfo_one_na <- wfo_na_check(
    wfo.result = wfo_one_nomatch$clean, 
    out.file = paste0(out.dir, "/wfo-one-na.csv"),
    verbose = verbose
  ) # removes NA scientificNames
  
  wfo_one_nomatch$clean = NULL # remove the clean to save memory
  
  wfo_one_dups <- wfo_duplication_check(
    wfo.result = wfo_one_na$clean, 
    out.file = paste0(out.dir, "/wfo-one-duplicates.csv"), 
    verbose = verbose
  ) # removes duplicate scientificNames
  
  wfo_one_na$clean = NULL # remove the clean to save memory
  
  wfo_one_res <- wfo_one_dups$clean
  wfo_one_dups$clean = NULL
  
  catn("Writing unique file to:", colcat(out_file, color = "output"))
  fwrite(wfo_one_res, out_file, bom = T)
  
  vebcat("WFO.one results recieved", color = "funSuccess")

  return(
    list(
      raw = wfo_one_checklist,
      clean = wfo_one_res,
      mismatch = wfo_one_mis$mismatch,
      nomatch = wfo_one_nomatch$nomatch,
      na = wfo_one_na$na,
      duplicate = wfo_one_dups$duplicate
    )
  )
}
