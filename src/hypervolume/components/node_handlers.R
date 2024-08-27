setup_node <- function(pro.dir, iteration, min.disk.space, verbose = FALSE) {
  current_disk_space <- get_disk_space("/export", units = "GB")

  # Make dir objects
  log_dir <- paste0(pro.dir, "/logs")
  init_log <- paste0(log_dir, "/init-log.txt")
  create_file_if(init_log)

  if (verbose) {
    try(init_log <- file(init_log, open = "at")) # for error handling before the process even starts
    sink(init_log, type = "output")
    sink(init_log, type = "message")
  }

  vebcat("Initating iteration", iteration, veb = verbose)
  vebcat("Creating directories.", veb = verbose)

  node_dir <- paste0(log_dir, "/nodes")
  locks_dir <- paste0(pro.dir, "/locks") # Recursive directory creation
  node_lock_dir <- paste0(locks_dir, "/node-iteration")

  # Create dirs of they do not exist
  create_dir_if(node_dir, node_lock_dir)

  vebcat("Creating files.", veb = verbose)

  # Make file objects
  warn_file <- paste0(log_dir, "/warning.txt")
  err_file <- paste0(log_dir, "/error.txt")
  node_it <- paste0(log_dir, "/node-iterations.txt")

  create_file_if(warn_file, err_file, node_it, keep = TRUE)

  while_msg <- FALSE

  while (TRUE) {
    if (!while_msg) {
      vebcat("The node is stuck in node-iterations traffic...", veb = verbose)
      while_msg <- TRUE
    }
    if (is.locked(node_lock_dir, lock.n = 1)) {
      Sys.sleep(3)
    } else {
      Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
      vebcat("Locking file.", veb = verbose)
      lock_node_it <- lock(node_lock_dir, lock.n = 1, paste0("Locked by iteration ", iteration))
      if (!is.null(lock_node_it) && file.exists(lock_node_it)) {
        break
      }
    }
  }
  while_msg <- FALSE

  vebcat("The node has escaped the node-iterations traffic.", veb = verbose)


  try(node_con <- file(node_it, open = "at"))
  identifier <- paste0("node", iteration)
  writeLines(identifier, node_con)
  close(node_con)

  # Check if the current disk space is less than the minimum required
  if (current_disk_space <= min.disk.space) stop("Insufficient disk space. Stopping processing.")

  catn("Init-log finished.")

  if (verbose) {
    sink(type = "message")
    sink(type = "output")
    close(init_log)
  }

  return(list(
    dir = node_dir,
    identifier = identifier,
    locks = node_lock_dir,
    lock.setup = lock_node_it,
    it.log = node_it,
    warn.file = warn_file,
    err.file = err_file
  ))
}

# --------------------------- Process node ---------------------------- #

process_node <- function(pro.dir, lock.dir, lock.setup, node.it.log, identifier, iteration, spec.list, columns.to.read, init.dt, hv.incl.threshold = 0.5, verbose, warn.file, err.file, fun) {
  log_dir <- paste0(pro.dir, "/logs")
  log_nodes <- paste0(log_dir, "/nodes")

  spec_filename <- spec.list[[iteration]]
  spec_name <- gsub(".csv", "", basename(spec_filename))
  spec_name <- trimws(spec_name)

  sp_node_log <- paste0(log_nodes, "/", iteration, "-", gsub(" ", config$species$file_separator, spec_name), "-log.txt")

  failed_its_log <- paste0(log_dir, "/failed-iterations.txt")
  highest_it_log <- paste0(log_dir, "/highest-iteration.txt")

  create_file_if(highest_it_log, failed_its_log, keep = TRUE)
  create_file_if(sp_node_log, failed_its_log)

  try(sp_log <- file(sp_node_log, open = "at"))
  sink(sp_log, type = "output")
  sink(sp_log, type = "message")

  spec <- fread(spec_filename, select = columns.to.read)

  vebcat("Unlocking file.", veb = verbose)

  unlock(lock.setup)

  spec <- spec[spec$occurrenceStatus == "PRESENT", ]

  spec <- spec[, "occurrenceStatus" := NULL]

  nobs <- nrow(spec)

  vebprint(head(spec, 2), verbose, text = "spec:")

  catn("Run Iteration", highcat(iteration))
  catn("Node Identifier:", identifier)
  catn("Using Species:", highcat(spec_name))
  catn("Using file:", highcat(spec_filename))
  catn("Species observations:", nobs)
  catn("Log observations:", log(nobs))
  catn()

  tryCatch(
    {
      # Condense data
      spec_condensed <- condense_taxons(spec.dt = spec)

      # cntry_condensed <- condense_country(spec.dt = spec)

      cntry_condensed <- find_wgsrpd_region(
        spec.dt = spec,
        projection = "longlat",
        longitude = "decimalLongitude",
        latitude = "decimalLatitude",
        wgsrpd.dir = "./resources/region/wgsrpd",
        wgsrpdlvl = "3",
        wgsrpdlvl.name = TRUE,
        unique = TRUE,
        verbose = verbose
      )

      # Subset data
      spec <- spec[, .(cleanName, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, countryCode, stateProvince, year)]
    },
    warning = function(w) warn(w, warn.file = warn.file, warn.txt = "Warning when condensing data", iteration = iteration),
    error = function(e) {
      close(sp_log)
      closeAllConnections()
      err(e, err.file = err.file, err.txt = "Error when condensing data", iteration = iteration)
    }
  )

  catn("initiating main function")

  tryCatch(
    {
      res <- fun(spec, spec_name)
    },
    error = function(e) {
      vebcat("Error occurred in main function for node", iteration, ":", e$message, color = "nonFatalError")
      sink(type = "message")
      sink(type = "output")
      close(sp_node_log)
    }
  )

  catn("Checking result and finishing up")

  cbmnd_res <- cbind(res, spec_condensed)

  # Duplicate the rows to match the length of cntry_condensed
  rep_cbmnd_res <- cbmnd_res[rep(seq_len(nrow(cbmnd_res)), nrow(cntry_condensed)), ]

  final_res <- cbind(rep_cbmnd_res, cntry_condensed)

  vebprint(nrow(final_res), verbose, "Final Result length:")
  vebprint(head(final_res, 2), verbose, "Final Result Head:")

  res <- check_hv_results(final_res, init.dt = init.dt, pro.dir, hv.incl.threshold = hv.incl.threshold, verbose = verbose)

  while_msg <- FALSE

  while (TRUE) {
    if (!while_msg) {
      catn("The node has entered node-iterations queue...")
      while_msg <- TRUE
    }
    if (is.locked(lock.dir, lock.n = 1)) {
      Sys.sleep(3)
    } else {
      Sys.sleep(runif(1, 0, 1)) # Add a random delay between 0 and 1 second
      lock_node_it <- lock(lock.dir, lock.n = 1, paste0("Locked by ", iteration, "_", gsub(" ", config$species$file_separator, spec_name)))
      if (!is.null(lock_node_it) && file.exists(lock_node_it)) {
        break
      }
    }
  }
  while_msg <- FALSE

  if (res) {
    fwrite(final_res, paste0(pro.dir, "/stats/stats.csv"), append = T, bom = T)
  } else {
    if (is.null(iteration) || is.na(iteration)) {
      catn("Iterations is NULL or NA:", iteration)
    } else {
      try(failed_con <- file(failed_its_log, open = "a"))
      writeLines(as.character(iteration), failed_con)
      close(failed_con)
    }

    try(err_con <- file(err.file, open = "a"))
    writeLines(paste0(
      "Hypervolume sequence failed the output check for node", iteration, " and species: ", gsub(" ", config$species$file_separator, spec_name)
    ), err_con)

    close(err_con)
  }
  # Get the node iterations and remove node from it
  try(node_con <- file(node.it.log, open = "r"))
  lines <- readLines(node_con)
  close(node_con)

  vebcat("Removing from node-iterations:", lines[which(grepl(paste0("^", identifier, "$"), lines))], veb = verbose)

  lines <- lines[-which(grepl(paste0("^", identifier, "$"), lines))] # Use regular expression to get exact match

  vebcat("Rewriting to node-iterations:", lines, veb = verbose)

  try(node_con <- file(node.it.log, open = "w"))

  writeLines(lines, node_con)

  close(node_con)

  # Get the highest iterations and write if higher
  try(high_con <- file(highest_it_log, open = "r"))
  highest_it <- as.integer(readLines(high_con))
  close(high_con)

  vebcat("Previous highest Iteration:", highest_it, veb = verbose)
  vebcat("This Iteration:", iteration, veb = verbose)

  if (length(highest_it) == 0 || highest_it < iteration) {
    vebcat("Setting new iteration to:", iteration)
    try(high_con <- file(highest_it_log, open = "w"))
    writeLines(as.character(iteration), high_con)
    close(high_con)

    if (verbose) {
      try(high_con <- file(highest_it_log, open = "r"))
      highest_it <- as.integer(readLines(high_con))
      close(high_con)

      vebcat("Actual highest iteration in file:", highest_it)
    }
  }

  catn("Unlocking node iterations.")

  unlock(lock_node_it)

  catn("Node finished successfully.")

  sink(type = "message")
  sink(type = "output")
  close(sp_log)
}
