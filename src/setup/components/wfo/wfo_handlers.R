#------------------------#
####     parallel     ####
#------------------------#

# Issue when parallel processing WFO.match -- oriSeq starts at 1 for each chunk
wfo_parallel <- function(checklist, cols, out.dir = "./", cores.max = 1, evals = NULL, counter = 1, verbose = FALSE) {
  catn("Analyzing: ", highcat(nrow(checklist)), "species using", highcat(cores.max),"cores.")
  if (!exists("system.speed.wfo", where = .GlobalEnv)) system.speed.wfo <- 3.5
  
  calculate_etc(
    timer.res = system.speed.wfo + 30, 
    cores = cores.max, 
    data.length = nrow(checklist),
    invisible = TRUE
  )
  
  catn("Sorting into chunks.")
  column <- cols$spec.name
  
  checklist[, speciesID := .GRP, by = column] # Create new OriSeqID for parallel
  
  setorder(checklist, speciesID)
  
  n_chunks <- ceiling(nrow(checklist) / cores.max)
  
  chunks <- split(checklist, rep(1:cores.max, each = n_chunks, length.out = nrow(checklist)))
  
  catn("Running WFO.match in", highcat(length(chunks)), "chunks with", highcat(n_chunks), "species in each chunk.")
  
  node_wfo_dir <- paste0(out.dir, "/wfo-match-nodes")
  create_dir_if(node_wfo_dir)
  
  catn("WFO.match progress can be found at:", colcat(node_wfo_dir, color = "output"))
  
  cl <- makeCluster(cores.max)
  
  vebcat("Including the necessary components in each core.", veb = verbose)
  
  cluster_params <- c(
    "chunks",
    "cols",
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
         WFO.file = WFO_file,
         spec.name = cols$spec.name,
         Genus = cols$Genus,
         Species = cols$Species,
         Infraspecific.rank = cols$Infraspecific.rank,
         Infraspecific = cols$Infraspecific,
         Authorship = cols$Authorship,
         verbose = verbose, 
         counter = counter
        )
       
       chunk_result <- data.table(chunk_result)
       
       chunk_result[, OriSeq := chunk$speciesID[match(chunk_result$OriSeq, seq_len(nrow(chunk)))]]
       chunk_result[, speciesID := NULL]
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
wfo_mismatch_check <- function(wfo.result, col.origin = "interimName", out.file = NULL, unchecked = FALSE, verbose = FALSE) {
  vebprint(unchecked, verbose, "Include unchecked results:")
  vebprint(names(wfo.result), verbose, "Input data table names:")
  vebprint(nrow(wfo.result), verbose, "Number of rows in the input data table:")
  
  dt <- copy(wfo.result)
  skip <- FALSE
  new_cols <- c("mismatch.old", "mismatch.scientific", "mismatch.any")
  if (all(new_cols %in% names(dt))) skip <- TRUE
  vebprint(skip, verbose, "Skip:")
  
  setnames(dt, old = col.origin, new = "input")
  
  # Build table subset
  cols_to_select <- c("verbatimName", "input")
  if (paste0(col.origin, "Authorship") %in% names(dt)) cols_to_select <- c(cols_to_select, paste0(col.origin, "Authorship"))
  cols_to_select <- c(cols_to_select, "scientificName", "genus.clean", "specificEpithet.clean", "other.clean", "fullName.clean", "extra.clean", "structure.clean",  "New.accepted", "Old.status", "Old.name", "name.clean", "Fuzzy.dist", "mismatch.old", "mismatch.scientific", "mismatch.any")
  if (unchecked) cols_to_select <- c(cols_to_select, "taxonomicStatus")
  
  vebprint(cols_to_select, verbose, "Columns to select:")
  vebprint(names(dt), verbose, "Updated data table names:")
  
  if (!skip) {
    dt[, scientificName := gsub("\\s+", " ", trimws(scientificName))] # Remove double spaces
    dt[, Old.name := gsub("\\s+", " ", trimws(Old.name))] # Remove double spaces
    dt[, input := { # Clean symbols & designations
      tmp <- clean_string(input, verbose)
      tmp <- clean_designations(tmp, config$species$standard_infraEpithets, verbose)
      clean_symbols(tmp, config$species$standard_symbols, verbose)
    }]
    
    dt[, c("genus.clean", "specificEpithet.clean", "name.clean", "other.clean", "fullName.clean", "extra.clean", "structure.clean") := {
      res <- lapply(seq_along(input), function(i) {
        cat("\rProcessing rows for mismatches:", i, "of", .N)
        flush.console()
        clean_spec_name(input[i], config$species$standard_symbols, config$species$standard_infraEpithets, verbose)
      });catn()
      list(
        vapply(res, function(x) x$genus, character(1)),
        vapply(res, function(x) x$specificEpithet, character(1)),
        vapply(res, function(x) x$cleanName, character(1)),
        vapply(res, function(x) x$other, character(1)),
        vapply(res, function(x) x$fullName, character(1)),
        vapply(res, function(x) x$extra, character(1)),
        vapply(res, function(x) x$structure, character(1))
      )
    }]
    
    # Remove hybrid, double spaces and trim
    non_spec <- c("MultiInfraspecificTaxon", "infraspecificTaxon")
    dt[, `:=`(
      name.clean = trimws(gsub(" +", " ", gsub("×", "", name.clean))),
      specificEpithet.clean = trimws(gsub(" +", " ", gsub("×", "", name.clean))),
      structure.clean = fifelse(structure.clean %in% non_spec, structure.clean, "species"),
      clean_name = gsub(" ", "", trimws(name.clean)),
      clean_old_name = gsub(" ", "", trimws(Old.name)),
      clean_scientific_name = gsub(" ", "", trimws(scientificName)),
      clean_species_name = gsub(" ", "", trimws(paste(genus.clean, specificEpithet.clean)))
    )]
    
    # If Old.name is not equal to name.clean/input without author names or as species name, then add to manual
    # if taxonRank is species, then build a species name
    # If taxonRank not a species, then use fullname
    dt <- dt[, `:=`(
      mismatch.old = New.accepted == TRUE & (
        clean_old_name != clean_name & (structure.clean == "species" & clean_old_name != clean_species_name)
      ),
      mismatch.scientific = New.accepted == FALSE & (
        clean_scientific_name != clean_name & (structure.clean == "species" & clean_scientific_name != clean_species_name)
      ),
      mismatch.any = NA
    )]
    
    # & Fuzzy.dist >= 5
    dt[, `:=` (
      mismatch.any = mismatch.old | mismatch.scientific,
      clean_name = NULL,
      clean_old_name = NULL,
      clean_scientific_name = NULL,
      clean_species_name = NULL
    )]
    
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
  
  # remove any identical input that are in the wfo output to not cause future issues when handling manually
  removed_species <- main_res[input %in% mis_res$input]
  mis_res <- rbind(mis_res, removed_species, fill = TRUE)
  main_res <- main_res[!input %in% mis_res$input]
  
  vebprint(mis_res, verbose, "Mismatching results:")
  catn("Found", highcat(nrow(mis_res)), "mismatching species.")
  
  if (nrow(mis_res) > 0) {
    if (!is.null(out.file)) {
      write_res <- unique(mis_res, by = "input")
      write_res <- write_res[, ..cols_to_select]
      setnames(write_res, old = "input", new = col.origin)
      catn("Writing to file:", colcat(out.file, color = "output"))
      fwrite(write_res, out.file, bom = TRUE)
    } 
  } 
  
  return(list(
    clean = setnames(main_res, old = "input", new = col.origin),
    mismatch = setnames(mis_res, old = "input", new = col.origin)
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

#------------------------#
####    WFO minimal   ####
#------------------------#

WFO.minimal <- function(WFO.file) {
  if (!inherits(try(file.info(WFO.file), silent = TRUE), "try-error")) {
    catn("Reading file")
    cols_to_read <- c("taxonID", "scientificName", "scientificNameAuthorship", "taxonRank", "nomenclaturalStatus", "taxonomicStatus")
    
    catn("Reading WFO data")
    wfo_data <- fread(WFO.file, select = cols_to_read)
    wfo_data <- wfo_data[ tolower(wfo_data$taxonomicStatus) != "synonym"]
    invisible(gc())
    
  } else {
    stop("Neeed to be a file input.")
  }
  
  return(wfo_data)
}

#------------------------#
####    WFO extract   ####
#------------------------#

WFO.extract <- function(x, WFO.data, verbose = FALSE) {
  spec_dt <- copy(x)
  
  spec_dt[, `:=` (
    scientificName.ORIG = scientificName, # Keep the original
    scientificName = NULL # Remove the column
  )]
  
  spec_dt[, `:=`(
    scientificName = fifelse(
      is.na(genus) | is.na(specificEpithet),
      NA_character_,
      fifelse(
        tolower(taxonRank) == "species",
        fifelse(
          is.na(scientificNameAuthorship) | scientificNameAuthorship == "",
          trimws(scientificName.ORIG),
          paste(trimws(scientificName.ORIG), trimws(scientificNameAuthorship))
        ),
        paste(trimws(genus), trimws(specificEpithet))
      )
    ),
    hasAuthorship = fifelse(
      is.na(genus) | is.na(specificEpithet),
      NA,
      fifelse(
        tolower(taxonRank) == "species" &
          !is.na(scientificNameAuthorship) &
          scientificNameAuthorship != "",
        TRUE,
        FALSE
      )
    )
  )]
  
  need_check <- spec_dt[tolower(taxonRank) != "species"]
  
  vebprint(need_check, verbose, "Need to be checked dt:")
  vebcat(highcat(nrow(spec_dt)), "scientificNames found")
  vebcat(highcat(nrow(need_check)), "non-species found and need to be checked")
  
  for (i in 1:nrow(spec_dt)) {
    cat("\rChecking species", i, "/", nrow(spec_dt))
    flush.console()
    
    tr <- spec_dt[i]$taxonRank
    if (tolower(tr) == "species") next
    
    row <- spec_dt[i]
    spec <- paste(trimws(row$genus), trimws(row$specificEpithet))
    spec_orig <- row$scientificName.ORIG
    spec_author <- row$scientificNameAuthorship
    spec_replace <- row$scientificName
    
    vebprint(paste0("^", spec), verbose, "spec:")
    vebprint(spec_orig, verbose, "spec_orig:")
    vebprint(spec_author, verbose, "spec_author:")
    vebprint(spec_replace, verbose, "spec_replace:")
    
    # Find all instances of the species
    checked <- WFO.data[grepl(paste0("^", spec), scientificName)]
    
    # Check instances with orig name USELLES
    name_index <- which(paste(checked$scientificName, checked$scientificNameAuthorship) %in% paste(spec_orig, spec_author))
    
    # Get length of species taxons
    spec_subset <- checked[which(checked$taxonRank == "species")]
    vebprint(spec_subset, verbose, "spec_subset")
    
    if (nrow(spec_subset) == 0) next
    
    # if more than one species
    if (nrow(spec_subset) > 1) {
      # save ids as characters
      ids <- spec_subset$taxonID
      
      # Find the smallest id
      smallest_id <- ids[which.min(as.numeric(gsub("wfo-", "", spec_subset$taxonID)))]
      vebprint(smallest_id, verbose, "smallest ID:")
      
      # Return new species name from smallest id
      spec_subset <- WFO.data[taxonID == smallest_id]
    }
    
    if (is.na(spec_subset$scientificNameAuthorship) || spec_subset$scientificNameAuthorship == "") {
      new_spec <- spec_subset$scientificName
      author <- FALSE
    } else {
      new_spec <- paste(spec_subset$scientificName, spec_subset$scientificNameAuthorship)
      author <- TRUE
    }
    
    vebprint(new_spec, verbose, "New species:")
    
    spec_dt[i]$scientificName <- new_spec
    spec_dt[i]$hasAuthorship <- author
  };catn()
  
  return(spec_dt)
}
