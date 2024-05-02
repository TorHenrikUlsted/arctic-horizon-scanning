check_syn_wfo <- function(checklist, column, folder, cores.max = 1, verbose = FALSE, counter = 1) {
  if (!"data.table" %in% class(checklist) && !"data.frame" %in% class(checklist)) {
    stop("The input data is not in the 'data.table' or 'data.frame' format.", print(class(checklist)))
  }
  
  vebcat("Initiating WFO synonym check", color = "funInit")
  
  cores.max <- min(nrow(checklist), cores.max)
  
  catn("Running the WFO synonym check with column:", colcat(column, color = "indicator"), "for table:")
  print(head(checklist, 3))
  
  catn("Analyzing: ", highcat(nrow(checklist)), "species using", highcat(cores.max),"cores.")
  
  if (!exists("system.speed.wfo", where = .GlobalEnv)) system.speed.wfo <- 3.5
  
  time.setup = 30
  
  eta <- (nrow(checklist) * system.speed.wfo) / cores.max + time.setup
  
  # Convert the estimated time to days, hours, minutes, and seconds
  days <- floor(eta / (24*60*60))
  hours <- floor((eta %% (24*60*60)) / (60*60))
  minutes <- floor((eta %% (60*60)) / 60)
  seconds <- round(eta %% 60, 2)
  
  
  vebcat(paste("Estimated wait time:", days, "days", hours, "hours", minutes, "minutes", round(seconds, 2), "seconds"), color = "timer")
  
  wfo_timer <- start_timer("wfo_match")
  
  catn("Sorting into chunks.")
  
  n_seq_chunk <- ceiling(nrow(checklist) / cores.max)
  
  # Create a list to store the row indices for each chunk - this is to edit the OriSeq in the end
  n_seq <- lapply(seq_len(cores.max), function(i) {
    ((i - 1) * n_seq_chunk + 1):min(i * n_seq_chunk, nrow(checklist))
  })
  
  chunks <- split(checklist, rep(1:cores.max, each = ceiling(nrow(checklist) / cores.max), length.out = nrow(checklist)))
  
  catn("Running WFO.match in", highcat(length(chunks)), "chunks with", highcat(ceiling(nrow(checklist) / cores.max)), "species in each chunk.")
  
  node_wfo_dir <- paste0(folder, "/wfo-match-nodes")
  create_dir_if(node_wfo_dir)
  
  catn("WFO.match progress can be found at:", colcat(node_wfo_dir, color = "indicator"))
  
  cl <- makeCluster(cores.max)
  
  export_vars <- c("node_wfo_dir")
  
  clusterExport(cl, export_vars, envir = environment())
  
  clusterEvalQ(cl, {
    library(parallel)
    library(WorldFlora)
    library(data.table)
    library(crayon)
    source("./src/utils/components/custom_colors.R")
    cc <- custom_colors()
    source("./src/utils/components/condition_handlers.R")
    source("./src/utils/components/time_tracker.R")
    source("./src/utils/components/file_managers.R")
    source("./src/utils/components/loader.R")
    WFO_file <- load_wfo()
  })
  
  
  wfo_checklist <- clusterApplyLB(cl, seq_along(chunks), function(i) {
    tryCatch({
      chunk <- chunks[[i]]
      
      log_file_out <- paste0(node_wfo_dir, "/", "node-", i, "-log.txt")
      create_file_if(log_file_out)
      
      try(log_file_out <- file(log_file_out, open = "at"))
      sink(log_file_out, type = "output")
      sink(log_file_out, type = "message")
      
      catn("Chunk number", i)
      catn("chunk length", nrow(chunk))
      
      matched_list <- WFO.match(spec.data = chunk, spec.name = column, WFO.file = WFO_file, verbose = verbose, counter = counter)
    }, error = function(e) {
      sink(type = "message")
      sink(type = "output")
      close(log_file_out)
      
      vebcat("Error when running WFO.match in iteration", i, "~ Stopping cluster and closing all connections.", color = "fatalError")
      end_timer(wfo_timer)
      stopCluster(cl)
      closeAllConnections()
      stop(e$message)
    })
    
    catn("Node finished.")
    
    sink(type = "message")
    sink(type = "output")
    close(log_file_out)
    
    invisible(gc())
    
    return(matched_list)
  })
  
  catn("Finishing up.")
  
  stopCluster(cl)
  
  wfo_checklist_bound <- rbindlist(wfo_checklist)
  
  wfo_checklist_bound <- set_df_utf8(wfo_checklist_bound)
  
  fwrite(wfo_checklist_bound, paste0(folder, "/wfo-match.csv"), bom = T)
  
  end_timer(wfo_timer)
  
  vebcat("WFO synonym check completed", color = "funSuccess")
  
  return(wfo_checklist)
}
