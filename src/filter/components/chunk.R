#------------------------#
####    Chunk data    ####
#------------------------#

get_file_stats <- function(filename, chunk.name, out.file) {
  if (!file.exists(out.file)) {
    catn("Calculating total rows and unique", chunk.name, "...")
    
    total_calc <- system_calc_uniq_and_rows(filename, chunk.name, sep = "\t")
    
    writeLines(as.character(total_calc), out.file)
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0(
        "GBIF download returned **", total_calc[[2]], "** unique ", chunk.name,
        " with **", total_calc[[1]], "** total rows"
      )
    )
  } else {
    total_calc <- as.integer(readLines(out.file))
  }
  
  return(total_calc)
}


# Main chunking function
chunk_data <- function(spec.occ, chunk.name = "species", chunk.column, chunk.dir, chunk.size = 1e7, iterations = NULL, verbose = TRUE) {
  is_file <- is.character(spec.occ) && file.exists(spec.occ)
  
  vebcat("Initiating", if(is_file) "file" else "loaded dataframe", "chunking protocol.", color = "funInit")
  
  chunk_path <- file.path(chunk.dir, chunk.name)
  create_dir_if(chunk_path)
  
  err_log <- file.path(chunk.dir, if(is_file) "file-error.txt" else "loaded-error.txt")
  create_file_if(err_log)
  
  if (is_file) {
    total_calc <- get_file_stats(spec.occ, chunk.name, file.path(chunk.dir, "total-calc.txt"))
    total_rows <- total_calc[[1]]
    unique_count <- total_calc[[2]]
    df_header <- fread(spec.occ, nrows = 0)
  } else {
    unique_count <- length(unique(spec.occ$species))
    total_rows <- nrow(spec.occ)
    df_header <- names(spec.occ)
  }
  
  catn("Can read", highcat(chunk.size), "chunks at a time")
  
  total_chunks <- ceiling(total_rows / chunk.size)
  
  if (is.null(iterations)) {
    iteration_file <- file.path(chunk.dir, if(is_file) "file-iteration.txt" else "loaded-iteration.txt")
    
    if (file.exists(iteration_file)) {
      i_start <- as.integer(readLines(iteration_file))
      if (length(i_start) == 0) i_start <- 0
      i_start <- i_start + 1 # Start on the next chunk instead of rerunning the last one
    } else {
      create_file_if(iteration_file)
      i_start <- 1
    }
    
    if (i_start > total_chunks) {
      return(vebcat("Data already chunked", color = "funSuccess"))
    }
    
    chunks <- seq(i_start, total_chunks)
  } else {
    chunks <- iterations
  }
  
  vebcat("The data has", highcat(unique_count), "unique", highcat(chunk.name), "and", highcat(total_rows), "total rows.")
  
  catn("Writing out files to:", colcat(chunk_path, color = "output"))
  
  catn("Running sequential chunking process\n")
  cat(sprintf("%7s | %14s | %18s | %13s | %17s\n",
              "Chunk", "Chunk n_rows", "Chunk unique rows", "Chunks total", "Chunks remaining"))
  
  for (i in chunks) {
    tryCatch({
      if (is_file) {
        data <- suppressWarnings(fread(spec.occ, skip = (i - 1) * chunk.size + 1, nrows = chunk.size, col.names = names(df_header), verbose = FALSE))
      }
      
      data <- data[species != ""]
      data[, (chunk.column) := species]
      
      data <- data[order(data[[chunk.column]]), ]
      species_list <- split(data, data[[chunk.column]])
      species_list <- species_list[!names(species_list) %in% c("", NA)]
      
      lapply(names(species_list), function(species) {
        filename <- gsub(" ", config$species$file_separator, species)
        filename <- sub(paste0(config$species$file_separator, "$"), "", filename)
        file_path <- file.path(chunk.dir, chunk.name, paste0(filename, ".csv"))
        
        fwrite(species_list[[species]], file_path, append = file.exists(file_path))
      })
      
      cat(
        sprintf(
          "\r%7d | %14d | %18d | %13d | %17d",
          i, nrow(data), uniqueN(data[[chunk.column]]), total_chunks, total_chunks - i
        )
      )
      
      flush.console()
      writeLines(as.character(i), file.path(chunk.dir, if(is_file) "file-iteration.txt" else "loaded-iteration.txt"))
    }, error = function(e) {
      try(err_con <- file(err_log, open = "at"))
      sink(err_con, type = "output")
      
      catn("\nError in iteration", i)
      print(e$message)
      
      sink(type = "output")
      close(err_con)
    }, finally = {
      closeAllConnections()
      invisible(gc())
    })
  };catn()
  
  vebcat("Chunking protocol completed successfully.", color = "funSuccess")
}

#------------------------#
####  rename chunks   ####
#------------------------#

identify_accepted <- function(dt, verbose = FALSE) {
  if (!any(tolower(dt$status) %in% "accepted" & 
           tolower(dt$rank) %in% "species")
  ) {
    
    vebcat(
      "Found", highcat(length(unique(dt$usageKey))), 
      "scientific names and", highcat(sum(tolower(dt$rank) == "species")),
      "species", veb = verbose
    )
    
    if (length(unique(dt$usageKey)) > 1) {
      # identify accepted usageKey
      accepted <- unique(dt$acceptedUsageKey)
      # remove all keys that are not accepted
      accepted <- accepted[!is.na(accepted)]
      
    } else {
      accepted <- unique(dt$usageKey)
    }
    
    # safety in case acceptedKey does not lead back to the accepted scientificName
    if (length(accepted) == 0) return(invisible())
    
  } else {
    vebcat("Found an accepted species", veb = verbose)
    return(invisible())
  }
  
  return(accepted)
}

process_accepted <- function(occ.data, checklist, accepted.keys, symbols, designations, verbose = FALSE) {
  
  if (length(accepted.keys) == 1) {
    tryCatch({
      # Get all potential names from checklist
      potential_names <- checklist[usageKey == accepted.keys | acceptedUsageKey == accepted.keys, ]
      accepted_name <- NULL
      
      # Try getting accepted name in order of preference
      if (nrow(potential_names[status == "ACCEPTED"]) > 0) {
        accepted_name <- potential_names[status == "ACCEPTED", scientificName][1]
      } else if (nrow(potential_names) > 0) {
        accepted_name <- potential_names[1, scientificName]
      }
      
      # Only try API if we couldn't get a name from checklist
      if (is.null(accepted_name) || is.na(accepted_name) || accepted_name == "") {
        dt <- tryCatch({
          gbif_retry(accepted.keys, "name_usage")
        }, error = function(e) {
          if (verbose) catn("GBIF API call failed:", e$message)
          return(NULL)
        })
        
        if (!is.null(dt) && !is.null(dt$data)) {
          accepted_name <- dt$data$scientificName
        }
      }
      
      # If we got a valid name, clean and return it
      if (!is.null(accepted_name) && !is.na(accepted_name) && accepted_name != "") {
        accepted_name <- clean_spec_name(accepted_name, symbols, designations)$cleanName
        return(list(accepted_name = accepted_name))
      }
      
      # If we got here, use original name from occurrence data
      original_name <- occ.data$scientificName[1]
      if (!is.na(original_name)) {
        original_name <- clean_spec_name(original_name, symbols, designations)$cleanName
        return(list(accepted_name = original_name))
      }
      
      if (verbose) catn("Could not resolve name for key", accepted.keys)
      return(NULL)
      
    }, error = function(e) {
      if (verbose) catn("Error processing key", accepted.keys, ":", e$message)
      return(NULL)
    })
  }
  
  result <- list()
  
  for (key in accepted.keys) {
    tryCatch({
      # Get all rows that could be relevant
      checklist_subset <- checklist[acceptedUsageKey == key | usageKey == key, ]
      
      if (nrow(checklist_subset) > 0) {
        accepted_name <- NULL
        
        # Try getting accepted name in order of preference
        if (nrow(checklist_subset[status == "ACCEPTED"]) > 0) {
          accepted_name <- checklist_subset[status == "ACCEPTED", scientificName][1]
        } else {
          accepted_name <- checklist_subset[1, scientificName]
        }
        
        # Try API only if necessary
        if (is.null(accepted_name) || is.na(accepted_name) || accepted_name == "") {
          dt <- tryCatch({
            gbif_retry(key, "name_usage")
          }, error = function(e) {
            if (verbose) catn("GBIF API call failed for key", key, ":", e$message)
            return(NULL)
          })
          
          if (!is.null(dt) && !is.null(dt$data)) {
            accepted_name <- dt$data$scientificName
          }
        }
        
        # If still no valid name, try using original name from data
        if (is.null(accepted_name) || is.na(accepted_name) || accepted_name == "") {
          original_rows <- occ.data[taxonKey == key, ]
          if (nrow(original_rows) > 0) {
            accepted_name <- original_rows$scientificName[1]
          }
        }
        
        if (!is.null(accepted_name) && !is.na(accepted_name) && accepted_name != "") {
          # Get all relevant keys
          relevant_keys <- unique(c(checklist_subset$usageKey, checklist_subset$acceptedUsageKey))
          relevant_keys <- relevant_keys[!is.na(relevant_keys)]
          
          # Get occurrence data
          occurrence_subset <- occ.data[taxonKey %in% relevant_keys, ]
          
          # Clean name
          accepted_name <- clean_spec_name(accepted_name, symbols, designations)$cleanName
          
          # Only process if we have data
          if (nrow(occurrence_subset) > 0) {
            occurrence_subset[, cleanName := accepted_name]
            result[[accepted_name]] <- occurrence_subset
            
            if (verbose) {
              catn("Processed", accepted_name, "with", 
                   nrow(checklist_subset), "taxa and", 
                   nrow(occurrence_subset), "occurrences")
            }
          }
        }
      }
    }, error = function(e) {
      if (verbose) {
        catn("Error processing key", key, ":", e$message)
      }
    })
  }
  
  # Return original data if we couldn't process anything
  if (length(result) == 0 && nrow(occ.data) > 0) {
    original_name <- occ.data$scientificName[1]
    if (!is.na(original_name)) {
      original_name <- clean_spec_name(original_name, symbols, designations)$cleanName
      occ.data[, cleanName := original_name]
      result[[original_name]] <- occ.data
      if (verbose) catn("Using original name", original_name, "for unresolved data")
    }
  }
  
  return(result)
}

handle_chunk_files <- function(processed.data, orig.file, file.sep = "_", verbose = FALSE) {
  if ("accepted_name" %in% names(processed.data)) {
    # Single key scenario
    filename <- file.path(dirname(orig.file), paste0(gsub(" ", file.sep, processed.data$accepted_name), ".csv"))
    
    file.rename(orig.file, filename)
    
  } else {
    # Multiple keys scenario
    new_files <- character()
    
    for (name in names(processed.data)) {
      filename <- file.path(dirname(orig.file), paste0(gsub(" ", file.sep, name), ".csv"))
      
      fwrite(processed.data[[name]], filename)
      
      if (file.exists(filename) && file.size(filename) > 0) {
        new_files <- c(new_files, filename)
        vebcat("Created file for", highcat(name), "with", 
                            highcat(nrow(processed.data[[name]])), "occurrences", veb = verbose)
      } else {
        vebcat("Failed to create file", highcat(filename), veb = verbose)
      }
    }
    
    if (length(new_files) == length(processed.data)) {
      file.remove(orig.file)
      vebcat("Removed original file", highcat(orig.file), veb = verbose)
    } else {
      vebcat("Not all new files were created successfully. Keeping original file.", veb = verbose)
    }
  }
}

rename_chunks <- function(spec.dir, symbols, designations, file.sep = "_", verbose = FALSE) {
  files <- list.files(spec.dir, full.names = TRUE)
  
  iteration_file <- file.path(dirname(spec.dir), "rename-iteration.txt")
  
  if (file.exists(iteration_file)) {
    i_start <- as.integer(readLines(iteration_file))
    if (length(i_start) == 0) i_start <- 0
    i_start <- i_start + 1 # Start on the next chunk instead of rerunning the last one
  } else {
    create_file_if(iteration_file)
    i_start <- 1
  }
  
  if (i_start > length(files)) return(catn("All files already renamed"))
  
  vebcat("Renaming files to accepted name", color = "funInit")
  
  for (i in i_start:length(files)) {
    cat("\rChecking file", highcat(i), "/", highcat(length(files)))
    
    # Check for name more than string length of 2, means already processed
    str_len <- length(str_split(basename(files[i]), file.sep)[[1]])
    
    if (str_len > 2) { # Means it has been processed before
      writeLines(as.character(i), iteration_file)
      next
    }
    
    dt <- fread(files[i])
    
    unique_names <- unique(dt$scientificName)
    
    checklist <- as.data.table(gbif_retry(unique_names, "name_backbone_checklist"))
    
    accepted <- identify_accepted(checklist)
    
    if (!is.null(accepted)) {
      processed_accepted <- process_accepted(dt, checklist, accepted, symbols, designations, verbose)
      
      handle_chunk_files(processed_accepted, files[i], file.sep, verbose)
    }
    
    writeLines(as.character(i), iteration_file)
    
  };catn()
  
  vebcat("All files renamed successfully", color = "funSuccess")
}

#------------------------#
####   Clean chunks   ####
#------------------------#

clean_chunks <- function(chunk.name, chunk.column, chunk.dir, sp_w_keys, iterations = NULL, verbose = F) {
  vebcat("Initiating chunk cleaning protocol.", color = "funInit")
  err_log <- paste0(chunk.dir, "/cleaner-error.txt")
  create_file_if(err_log)
  
  # File to store the last iteration number
  last_iteration_file <- paste0(chunk.dir, "/cleaner-iteration.txt")
  missing_sp_file <- paste0(chunk.dir, "/missing-species.txt")
  
  # Get the list of all chunk files
  chunk_files <- list.files(paste0(chunk.dir, "/", chunk.name), pattern = "*.csv", full.names = TRUE)
  
  # Total number of chunk files
  n_total <- length(chunk_files)
  
  if (file.exists(last_iteration_file)) {
    if (is.null(iterations)) {
      iterations <- as.integer(readLines(last_iteration_file))
      
      if (iterations >= n_total) {
        return(vebcat("This data has already been cleaned.", color = "funSuccess"))
      }
    }
  }
  
  
  if (is.null(iterations) || is.na(iterations)) {
    iterations <- seq_along(chunk_files)
  }
  
  vebprint(chunk_files, text = "chunk_files:", veb = verbose)
  
  j <- 0
  
  catn("Cleaning chunk files")
  cat(sprintf("%6s | %10s | %10s \n", "total", "iteration", "Remaining"))
  
  for (i in iterations) {
    cat(sprintf("\r%6.0f | %10.0f | %10.0f", n_total, i, n_total - i))
    flush.console()
    
    strings <- gsub("\\s+", "", sp_w_keys$species)
    
    file_name <- chunk_files[i]
    
    name <- basename(file_name)
    
    name <- gsub(".csv", "", name)
    
    name <- gsub(config$species$file_separator, " ", name)
    
    name_nospace <- gsub("\\s+", "", name)
    
    # Fuzzy match 1 letter -- birngs too many errors
    # matches <- agrepl(name_nospace, strings, max.distance = 0.001)
    
    if (name_nospace %in% strings) {
      next
    } else {
      vebcat("\nRemoving", highcat(name), veb = verbose)
      file.remove(file_name)
    }
    
    j + 1
  }
  catn()
  
  try(it_con <- open(last_iteration_file, "w"))
  writeLines(as.character(j), it_con)
  close(it_con)
  
  chunk_files <- gsub("\\s+", "", gsub(config$species$file_separator, "", gsub(".csv", "", basename(list.files(paste0(chunk.dir, "/", chunk.name), pattern = "*.csv", full.names = TRUE)))))
  
  missing_sp <- chunk_files[!(chunk_files %in% gsub("\\s+", "", sp_w_keys$species))]
  
  if (length(missing_sp) > 0) {
    vebcat("These species are not found as exact matches in the species keys", color = "indicator")
    try(sp_con <- open(missing_sp_file, "w"))
    writeLines(missing_sp, sp_con)
    close(sp_con)
  } else {
    catn("All chunk files are present in sp_w_keys$species.")
  }
  
  list_md <- list.files(paste0(chunk.dir, "/", chunk.name), pattern = "*.csv", )
  
  mdwrite(
    config$files$post_seq_md,
    text = paste0(
      "Number of species after cleaning: **", length(list_md), "**"
    )
  )
  
  vebcat("Chunk cleaning protocol completed successfully.", color = "funSuccess")
}
