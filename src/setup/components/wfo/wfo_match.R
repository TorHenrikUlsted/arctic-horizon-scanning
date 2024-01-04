check_syn_wfo <- function(checklist, column, folder, max.cores, verbose, counter) {
  if (!"data.table" %in% class(checklist) && !"data.frame" %in% class(checklist)) {
    stop("The input data is not in the 'data.table' or 'data.frame' format.", print(class(checklist)))
  }
  
  max.cores <- min(nrow(checklist), max.cores)
  
  cat("Running the WFO synonym check with column:", cc$aquamarine(column), "for table: \n")
  print(head(checklist, 3))
  
  cat("Analyzing: ", magenta(nrow(checklist)), "species using", magenta(max.cores),"cores. \n")
  
  if (!exists("time.const", where = .GlobalEnv)) {
    time.const <- 3.5
    time.setup = 30
  }
  eta <- (nrow(checklist) * time.const) / max.cores + time.setup
  
  # Convert the estimated time to days, hours, minutes, and seconds
  days <- floor(eta / (24*60*60))
  hours <- floor((eta %% (24*60*60)) / (60*60))
  minutes <- floor((eta %% (60*60)) / 60)
  seconds <- round(eta %% 60, 2)
  
  
  cat(magenta(paste("Estimated wait time:", days, "days", hours, "hours", minutes, "minutes", round(seconds, 2), "seconds\n")))

  wfo_timer <- start_timer("wfo_match")
  
  cat("Sorting into chunks. \n")
  
  n_seq_chunk <- ceiling(nrow(checklist) / max.cores)
  
  # Create a list to store the row indices for each chunk - this is to edit the OriSeq in the end
  n_seq <- lapply(seq_len(max.cores), function(i) {
    ((i - 1) * n_seq_chunk + 1):min(i * n_seq_chunk, nrow(checklist))
  })
  
  chunks <- split(checklist, rep(1:max.cores, each = ceiling(nrow(checklist) / max.cores), length.out = nrow(checklist)))
  
  cat("Running WFO.match in", cc$lightSteelBlue(length(chunks)), "chunks with", cc$lightSteelBlue(ceiling(nrow(checklist) / max.cores)), "species in each chunk. \n")
  
  node_hv_dir <- paste0(folder, "/wfo-match-nodes")
  create_dir_if(node_hv_dir)
  
  cat("WFO.match progress can be found at:", yellow(node_hv_dir), "\n")
  
  cl <- makeCluster(max.cores)
  
  clusterEvalQ(cl, {
    library(WorldFlora)
    library(data.table)
    source("./src/utils/components/get_wfo_backbone.R")
  })
  
  
  wfo_checklist <- clusterApplyLB(cl, seq_along(chunks), function(i) {
    chunk <- chunks[[i]]
    
    log_file_out <- paste0(node_hv_dir, "/", "node-", i, "-log.txt")
    
    if (!file.exists(log_file_out)) {
      file.create(log_file_out)
    } else {
      file.remove(log_file_out)
      file.create(log_file_out)
    }
    
    try(log_file_out <- file(log_file_out, open = "at"))
    sink(log_file_out, type = "output")
    try(log_file_out <- file(log_file_out, open = "at"))
    sink(log_file_out, type = "message")
    
    cat("Chunk number", i, "\n")
    cat("chunk length", nrow(chunk), "\n")
    
    tryCatch({
      matched_list <- WFO.match(spec.data = chunk, spec.name = column, WFO.file = WFO_file, verbose = verbose, counter = counter)
    }, error = function(e) {
      cat("Error when running WFO.match in iteration", i, "\n")
      cat(e)
    })
    
    cat("Node finished. \n")
    
    sink(type = "message")
    sink(type = "output")
    close(log_file_out)
    
    invisible(gc())
    
    return(matched_list)
  })
  
  cat("Finishing up \n")
  
  stopCluster(cl)
  
  wfo_checklist_bound <- rbindlist(wfo_checklist)
  
  wfo_checklist_bound <- set_df_utf8(wfo_checklist_bound)

  fwrite(wfo_checklist_bound, paste0(folder, "/wfo-match.csv"), bom = T)

  end_timer(wfo_timer)

  cat(cc$aquamarine("WFO synonym check completed \n"))
  
  return(wfo_checklist)
}
