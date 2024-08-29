#------------------------#
####     parallel     ####
#------------------------#

# Issue when parallel processing WFO.match -- oriSeq starts at 1 for each chunk
wfo_parallel <- function(checklist, column, out.dir = "./", cores.max = 1, evals = NULL, counter = 1, verbose = FALSE) {
  catn("Analyzing: ", highcat(nrow(checklist)), "species using", highcat(cores.max),"cores.")
  if (!exists("system.speed.wfo", where = .GlobalEnv)) system.speed.wfo <- 3.5
  
  calculate_etc(
    timer.res = system.speed.wfo + 30, 
    cores = cores.max, 
    data.length = nrow(checklist),
    invisible = TRUE
  )
  
  catn("Sorting into chunks.")
  
  checklist[, speciesID := .GRP, by = column]  # column is your species name column
  
  setorder(checklist, speciesID)
  
  n_chunks <- ceiling(nrow(checklist) / cores.max)
  
  chunks <- split(checklist, rep(1:cores.max, each = n_chunks, length.out = nrow(checklist)))
  
  catn("Running WFO.match in", highcat(length(chunks)), "chunks with", highcat(n_chunks), "species in each chunk.")
  
  node_wfo_dir <- paste0(out.dir, "/wfo-match-nodes")
  create_dir_if(node_wfo_dir)
  
  catn("WFO.match progress can be found at:", colcat(node_wfo_dir, color = "indicator"))
  
  cl <- makeCluster(cores.max)
  
  vebcat("Including the necessary components in each core.", veb = verbose)
  
  cluster_params <- c(
    "chunks",
    "column",
    "node_wfo_dir",
    "cores.max",
    "evals",
    "verbose"
  )
  
  clusterExport(cl, cluster_params, envir = environment())
  
  clusterEvalQ(cl, {
    lapply(evals$packages, require, character.only = TRUE)
    lapply(evals$source, function(file) {
      tryCatch({
        source(file)
      }, error = function(e) {
        warning(paste("Error sourcing file:", file, "\nError message:", e$message))
      })
    })
  })
  
 tryCatch({
   wfo_parallel_result <- clusterApplyLB(cl, seq_along(chunks), function(i) {
     tryCatch({
       chunk <- chunks[[i]]
       
       log_file_out <- paste0(node_wfo_dir, "/", "node-", i, "-log.txt")
       create_file_if(log_file_out)
       
       try(log_file_out <- file(log_file_out, open = "at"))
       sink(log_file_out, type = "output")
       sink(log_file_out, type = "message")
       
       catn("Chunk number", i)
       catn("chunk length", nrow(chunk))
       catn("Is data.table:", is.data.table(chunk))
       catn("Is data.frame:", is.data.frame(chunk))
       catn("Loading WFO backbone")
       WFO_file <- load_wfo()
       
       chunk_result <- WFO.match(
         spec.data = chunk, 
         spec.name = column, 
         WFO.file = WFO_file, 
         verbose = verbose, 
         counter = counter
        )
       
       chunk_result <- data.table(chunk_result)
       
       chunk_result[, OriSeq := chunk$speciesID[match(chunk_result$OriSeq, seq_len(nrow(chunk)))]]
     }, error = function(e) {
       message("Error when running WFO.match in iteration ", i, " ~ Stopping cluster and closing all connections.")
       stop(e)
     }, finally = {
       catn("Cleaning up node.")
       sink(type = "message")
       sink(type = "output")
       close(log_file_out)
     })
     catn("Node finished.")
     return(chunk_result)
   })
 }, error = function(e) {
   vebcat("Error when running WFO.match in iteration ", i, " ~ Stopping cluster and closing all connections.", color = "fatalError")
   stop(e$message)
 }, finally = {
   catn("Cleaning up connections")
   stopCluster(cl)
   closeAllConnections()
   invisible(gc())
 })
  
  wfo_result <- rbindlist(wfo_parallel_result)
  
  wfo_result[, Subseq := seq_len(.N), by = OriSeq]
  
  return(wfo_result)
}

#------------------------#
####   mismatch case  ####
#------------------------#
wfo_mismatch_check <- function(wfo.result, col.origin = "rawName", out.file = NULL, unchecked = FALSE, verbose = FALSE) {
  vebprint(unchecked, verbose, "Include unchecked results:")
  vebprint(names(wfo.result), verbose, "Input data table names:")
  vebprint(nrow(wfo.result), verbose, "Number of rows in the input data table:")
  
  dt <- copy(wfo.result)
  skip <- FALSE
  new_cols <- c("mismatch.input", "mismatch.old", "mismatch.scientific", "mismatch.any")
  if (all(new_cols %in% names(dt))) skip <- TRUE
  vebprint(skip, verbose, "Skip:")
  
  input_orig <- paste0(col.origin, ".ORIG")
  setnames(dt, old = c(col.origin, input_orig), new = c("input", "input.ORIG"))
  
  cols_to_select <- c("input", "input.ORIG", "scientificName", "New.accepted", "Old.status", "Old.name", "Fuzzy.dist", "mismatch.input", "mismatch.old", "mismatch.scientific", "mismatch.any")
  if (unchecked) cols_to_select <- c(cols_to_select, "taxonomicStatus")
  
  vebprint(cols_to_select, verbose, "Columns to select:")
  vebprint(names(dt), verbose, "Updated data table names:")
  
  if (!skip) {
    dt[, scientificName := gsub("\\s+", " ", trimws(scientificName))] # Remove double spaces
    dt[, input.ORIG := gsub("\\s+", " ", trimws(input.ORIG))] # Remove double spaces
    dt[, input := gsub("\\s+", " ", trimws(input))] # Remove double spaces
    dt[, Old.name := gsub("\\s+", " ", trimws(Old.name))] # Remove double spaces
    
    dt <- dt[, `:=`(
      mismatch.input = input != input.ORIG,
      mismatch.old = New.accepted == TRUE & input.ORIG != Old.name,
      mismatch.scientific = New.accepted == FALSE & input.ORIG != trimws(scientificName),
      mismatch.any = NA
    )]
    
    dt[, mismatch.any := mismatch.input | mismatch.old | mismatch.scientific]
  } else {
    catn("Using preexisting columns")
  }
  
  if (unchecked) {
    main_res <- dt[(dt$mismatch.any == FALSE | tolower(dt$taxonomicStatus) != "unchecked")]
    mis_res <- dt[(dt$mismatch.any == TRUE | tolower(dt$taxonomicStatus) == "unchecked")]
  } else {
    main_res <- dt[(dt$mismatch.any == FALSE)]
    mis_res <- dt[(dt$mismatch.any == TRUE)]
  }
  # remove any identical input.ORIG that are in the wfo output to not cause future issues when handling manually
  removed_species <- main_res[input.ORIG %in% mis_res$input.ORIG]
  mis_res <- rbind(mis_res, removed_species, fill = TRUE)
  main_res <- main_res[!input.ORIG %in% mis_res$input.ORIG]
  
  vebprint(mis_res, verbose, "Mismatching results:")
  catn("Found", highcat(nrow(mis_res)), "mismatching species.")
  
  if (nrow(mis_res) > 0) {
    if (!is.null(out.file)) {
      write_res <- unique(mis_res, by = "input.ORIG")
      write_res <- write_res[, ..cols_to_select]
      setnames(write_res, old = c("input", "input.ORIG"), new = c(col.origin, input_orig))
      catn("Writing to file:", colcat(out.file, color = "output"))
      fwrite(write_res, out.file, bom = TRUE)
    } 
  } 
  
  return(list(
    clean = setnames(main_res, old = c("input", "input.ORIG"), new = c(col.origin, input_orig)),
    mismatch = setnames(mis_res, old = c("input", "input.ORIG"), new = c(col.origin, input_orig))
  ))
}

#------------------------#
####   nomatch case   ####
#------------------------#
wfo_nomatch_check <- function(wfo.result, out.file = NULL, verbose = FALSE) {
  res <- copy(wfo.result)
  
  nomatch <- res[tolower(res$One.Reason) == "no match found"]
  res <- res[tolower(res$One.Reason) != "no match found"]
  
  vebprint(nrow(nomatch), verbose, "Number of nomatches")
  vebprint(nrow(res), verbose, "Number of rows in result")
  
  if (nrow(nomatch) > 0 && !is.null(out.file)) {
    vebcat("Missing matches found for", highcat(nrow(nomatch)), "species, manually check these.", color = "nonFatalError")
    
    catn("Writing out to:", colcat(out.file, color = "output"))
    
    fwrite(nomatch, out.file, bom = TRUE)
  }
  
  catn("Removed", highcat(nrow(nomatch)), "missing matches from WFO result.")
  
  return(list(
    clean = res,
    nomatch = nomatch
  ))
}

#------------------------#
####  duplicate case  ####
#------------------------#
wfo_duplication_check <- function(wfo.result, out.file = NULL, verbose = FALSE) {
  res <- copy(wfo.result)
  
  dups <- res[duplicated(res$scientificName)]
  res <- res[!duplicated(res$scientificName)]
  
  vebprint(dups, verbose, "Duplications dt:")
  vebprint(wfo.result, verbose, "Result dt:")
  
  if (nrow(dups) > 0 && !is.null(out.file)) {
    catn("Found", highcat(nrow(dups)), "duplicated species from the WFO.one result.")
    
    catn("Writing out to:", colcat(out.file, color = "output"))
    
    fwrite(dups, out.file, bom = T)
  }
  
  catn("Removed", highcat(nrow(dups)), "species from the WFO result.")
  
  return(list(
    clean = res,
    duplicate = dups
  ))
}

#------------------------#
####      NA case     ####
#------------------------#
wfo_na_check <- function(wfo.result, out.file = NULL, verbose = FALSE) {
  res <- copy(wfo.result)
  
  res_na <- res[is.na(res$scientificName)]
  res <- res[!is.na(res$scientificName)]
  
  catn("Found", highcat(nrow(res_na)), "NA scientificNames")
  
  if (nrow(res_na) > 0) {
    if (!is.null(out.file)) {
      catn("Writing to file:", colcat(out.file, color = "output"))
      fwrite(res, out.file, bom = T)
    }
  }
  
  return(list(
    clean = res,
    na = res_na
  ))
}

