setup_node <- function(pro.dir, iteration, min.disk.space, verbose = FALSE) {
  current_disk_space <- get_disk_space("/export", units = "GB")
  
  # Make dir objects
  log_dir <- paste0(pro.dir, "/logs")
  init_log <- paste0(log_dir, "/init-log.txt")
  create_file_if(init_log)
  
  node <- paste0("Node-", iteration)
  
  try(init_log <- file(init_log, open = "at")) # for error handling before the process even starts
  sink(init_log, type = "output")
  sink(init_log, type = "message")
  
  node_timer <- start_timer(node)
  
  vebcat("Creating directories.", veb = verbose)
  
  node_dir <- paste0(log_dir, "/nodes")
  locks_dir <- paste0(pro.dir, "/locks") # Recursive directory creation
  node_lock_dir <- paste0(locks_dir, "/node-iteration")
  
  # Create dirs of they do not exist
  create_dir_if(c(node_dir, node_lock_dir))
  
  vebcat("Creating files.", veb = verbose)
  
  # Make file objects
  warn_file <- paste0(log_dir, "/warning.txt")
  err_file <- paste0(log_dir, "/error.txt")
  node_it <- paste0(log_dir, "/node-iterations.txt")
  while_msg <- FALSE
  
  vebcat("Locking file.", veb = verbose)
  
  while (TRUE) {
    if (!while_msg) {
      vebcat("The node is stuck in node-iterations traffic...", veb = verbose)
      while_msg <- TRUE
    }
    if (is.locked(lock_node_dir, lock.n = 1)) {
      Sys.sleep(1)
    } else {
      Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
      lock_node_it <- lock(lock_node_dir, lock.n = 1)
      if (!is.null(lock_node_it) && file.exists(lock_node_it)) {
        break
      }
    }
  }
  while_msg <- FALSE
  
  vebcat("Checking for locked file.", veb = verbose)
  
  node_con <- file(node_it, open = "at")
  identifier <- paste0("node", iteration)
  writeLines(identifier, node_con)
  close(node_con)
  
  unlock(lock_node_it)
  
  # Check if the current disk space is less than the minimum required
  if (current_disk_space <= min.disk.space) stop("Insufficient disk space. Stopping processing.")
  
  catn("Init-log finished.")
  
  sink(type = "message")
  sink(type = "output")
  close(init_log)
  
  return(list(
    identifier = identifier,
    dir = node_dir,
    locks = node_lock_dir, 
    warn = warn_file,
    err = err_file
  ))
}

# --------------------------- Process node ---------------------------- #

process_node <- function(pro.dir, iteration, identifier, spec.list, columns.to.read, init.dt, hv.incl.threshold = 0.5, verbose, warn, err, fun) {
  log_dir <- paste0(pro.dir, "/logs")
  
  spec_filename <- spec.list[[iteration]]
  spec_name <- gsub(".csv", "", basename(spec_filename))
  spec_name <- trimws(spec_name)
  
  sp_node_log <- paste0(log_dir, "/", iteration, "_", gsub(" ", "-", spec_name), "-log.txt")
  create_file_if(sp_node_log)
  
  try(sp_node_log <- file(sp_node_log, open = "at"))
  sink(sp_node_log, type = "output")
  sink(sp_node_log, type = "message")
  
  spec <- fread(spec_filename, select = columns.to.read)
  
  spec <- spec[spec$occurrenceStatus == "PRESENT", ]
  
  spec <- spec[, "occurrenceStatus" := NULL]
  
  nobs <- nrow(spec)
  
  vebprint(head(spec, 2), text = "spec:")
  
  catn("Run Iteration", highcat(iteration))
  catn("Node Identifier:", identifier, "\n")
  catn("Using Species:", highcat(spec.name))
  catn("Using Species:", highcat(spec_filename))
  catn("Species observations:", nobs)
  
  # Condense data
  spec_condensed <- condense_taxons(spec.dt = spec)
  
  cntry_condensed <- condense_country(spec.dt = spec)
  
  tryCatch({
    res <- fun(spec, spec_name)
  },
  error = function(e) {
    vebcat("Error occurred in node", iteration, ":", e$message, color = "nonFatalError")
    sink(type = "message")
    sink(type = "output")
    close(sp_node_log)
  })
  
  cbmnd_res <- cbind(res, spec_condensed)
  
  # Duplicate the rows to match the length of cntry_condensed
  rep_cbmnd_res <- cbmnd_res[rep(seq_len(nrow(cbmnd_res)), nrow(cntry_condensed)), ]
  
  final_res <- cbind(rep_cbmnd_res, cntry_condensed)
  
  vebprint(final_res, verbose, "Final Result:")
  
  res <- check_results(final_res, pro.dir, init.dt = init.dt, verbose = verbose)
  
  while (TRUE) {
    if (!while_msg) {
      catn("The node has entered node-iterations queue...")
      while_msg <- TRUE
    }
    if (is.locked(lock_node_dir, lock.n = 1)) {
      Sys.sleep(1)
    } else {
      Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
      lock_node_it <- lock(lock_node_dir, lock.n = 1)
      if (!is.null(lock_node_it) && file.exists(lock_node_it)) {
        break
      }
    }
  }
  while_msg <- FALSE
  
  if (res) {
    fwrite(final_res, paste0(pro.dir, "/stats.csv"), append = T, bom = T)
  } else {
    err_con <- file(err, open = "a")
    
    writeLines(paste0(
      "Hypervolume sequence failed the output check for node", iteration
    ), err_con)
    
    close(err_con)
  }
  
  node_con <- file(node_it, open = "r")
  
  lines <- readLines(node_con)
  
  close(node_con)
  
  vebcat("Removing from node-iterations:", lines[which(grepl(paste0("^", identifier, "$"), lines))], veb = verbose)
  
  lines <- lines[-which(grepl(paste0("^", identifier, "$"), lines))] # Use regular expression to get exact match
  
  vebcat("Rewriting to node-iterations:", lines, veb = verbose)
  
  node_con <- file(node_it, open = "w")
  
  writeLines(lines, node_con)
  
  close(node_con)
  
  catn("Unlocking node iterations.")
  
  unlock(lock_node_it)
  
  catn("Node finished successfully.")
  
  end_timer(node_timer)
  
  sink(type = "message")
  sink(type = "output")
  close(sp_node_log)
}