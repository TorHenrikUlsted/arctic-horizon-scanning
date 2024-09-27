#####################
# Chunk from file
#####################

chunk_file <- function(file.path, approach = "precautionary", cores.max = 1, chunk.name, chunk.column, chunk.dir, chunk.size = 1e6, iterations = NULL, verbose = T) {
  vebcat("Initiating file chunking protocol.", color = "funInit")

  if (!file.exists(file.path)) {
    vebcat("Cannot find file.path, reading", file.path, color = "nonFatalError")
  } else {
    catn("Using file", highcat(file.path))
  }
  
  chunk_path <- paste0(chunk.dir, "/", chunk.name)
  create_dir_if(chunk_path)
  
  calc_file <- paste0(chunk.dir, "/total-calc.txt")
  err_log <- paste0(chunk.dir, "/file-error.txt")
  create_file_if(err_log)

  # Calculate the total number of chunks
  if (!file.exists(calc_file)) {
    catn("Calculating total rows and unique", chunk.name, "...")
    total_calc <- system_calc_uniq_and_rows(file.path, chunk.name, sep = "\t")
    writeLines(as.character(total_calc), calc_file)
    mdwrite(
      config$files$post_seq_md,
      text = paste0(
        "GBIF download returned **", total_calc[[2]], "** unique ", chunk.name,
        " with **", total_calc[[1]], "** total rows"
      )
    )
  } else {
    total_calc <- as.integer(readLines(calc_file))
  }
  
  total_rows <- total_calc[[1]]
  unique_count <- total_calc[[2]]
  
  if (is.null(iterations)) {
    
    vebcat("The file has", highcat(unique_count), "unique", highcat(chunk.name), "and", highcat(total_rows), "total rows.")
    vebcat("chunk.size: ", highcat(chunk.size), veb = verbose)
    
    total_chunks <- ceiling(total_rows / chunk.size)
        
    vebcat("Total chunks:", highcat(total_chunks), veb = verbose)

    iteration_file <- paste0(chunk.dir, "/file-iteration.txt")

    if (file.exists(iteration_file)) {
      i_start <- as.integer(readLines(iteration_file))
      if (length(i_start) == 0) i_start <- 0
      i_start <- i_start + 1 # Start on the next chunk instead of rerunning the last one
    } else {
      create_file_if(iteration_file)
      i_start <- 1
    }
  }
  
  catn("Initiating chunking protocol from", highcat(i_start), "/", highcat(total_chunks), "total chunks")

  if (i_start > total_chunks) {
    return(vebcat("Data already chunked", color = "funSuccess"))
  }

  if (is.vector(chunk.column) && length(chunk.column) > 1) {
    vebcat("Combining columns:", highcat(chunk.column), veb = verbose)
  }
  
  # Create progress directory
  progress_dir <- file.path(chunk.dir, "nodes")
  create_dir_if(progress_dir)
  
  progress_file <- file.path(progress_dir, paste0("node_", seq_len(cores.max), ".txt"))
  file.create(progress_file)
  
  df_header <- fread(file.path, nrows = 0)
  vebprint(df_header, verbose, "data table header:")

  catn("Writing out species files to:", colcat(chunk_path, color = "output"))

  process_chunk <- function(chunk_info) {
    i <- chunk_info$chunk_id
    core_id <- chunk_info$core_id
    
    if (cores.max > 1) {
      con <- file(paste0(progress_dir, "/node_", core_id, ".txt"), "w")
      sink(con)
      catn(paste(i, "/", total_chunks))
    }
    
    tryCatch({
      # Read chunk with suppressed because of incorrect quotation marks in some cases
      data <- suppressWarnings(fread(file.path, skip = (i - 1) * chunk.size + 1, nrows = chunk.size, col.names = names(df_header), verbose = FALSE))
      
      # Process data
      data <- select_species_approach(
        dt = data,
        approach = approach,
        col.name = chunk.column,
        custom.list = config$species$taxonRank_infraEpithet,
        verbose = verbose
      )
      
      # Split data by species
      data <- data[order(data[[chunk.column]]), ]
      species_list <- split(data, data[[chunk.column]])
      species_list <- species_list[!names(species_list) %in% c("", NA)]
      
      # Write each species data to file
      lapply(names(species_list), function(species) {
        filename <- gsub(" ", config$species$file_separator, species)
        filename <- sub(paste0(config$species$file_separator, "$"), "", filename)
        file_path <- file.path(chunk.dir, chunk.name, paste0(filename, ".csv"))
        
        fwrite(species_list[[species]], file_path, append = file.exists(file_path))
      })
      
    if (cores.max == 1) {
      cat(
        sprintf(
          "\r%7.0f | %14.0f | %18.0f | %13s | %13.0f | %17.0f",
          i, nrow(data), length(unique(data[[chunk.column]])), class(species_list), total_chunks, ceiling((total_rows - (i * chunk.size)) / chunk.size)
        )
      )
      flush.console()
      writeLines(as.character(i), iteration_file)
    }
    }, error = function(e) {
      err_con <- try(err_log, )
      try(err_con <- file(err_log, open = "at"))
      sink(err_con, type = "output")
      
      catn("\nError in iteration", i)
      print(e$message)
      
      sink(type = "output")
      close(err_con)
    }, finally = {
      if (cores.max > 1) {
        catn("Cleaning up connections")
        sink(output)
        
      }
      closeAllConnections()
      invisible(gc())
    })
  }
  
  # Prepare chunk information
  chunks <- seq(i_start, total_chunks)
  chunk_info <- lapply(seq_along(chunks), function(i) list(chunk_id = chunks[i], core_id = (i - 1) %% cores.max + 1))
  
  if (cores.max > 1) {
    catn("Setting up cores")
    # Process chunks in parallel
    cl <- makeCluster(cores.max)
    on.exit(stopCluster(cl))
    
    # Export all necessary functions and variables to the cluster
    clusterExport(cl, c("select_species_approach", "config", "progress_dir", 
                        "total_chunks", "chunk.dir", "chunk.name", "file.path", "chunk.size", 
                        "names", "df_header", "approach", "chunk.column"), envir = environment())
    clusterEvalQ(cl, {
      source("./src/utils/utils.R")
      load_utils(parallel = TRUE)
    })
    
    catn("Running parallel chunking process")
    catn("Writing out process to:", colcat(progress_dir, color = "output"))
    results <- parLapply(cl, chunk_info, process_chunk)
  } else {
    catn("Running sequential chunking process\n")
    cat(
      sprintf(
        "%7s | %14s | %18s |  %12s | %13s | %17s\n",
        "Chunk", "Chunk n_rows", "Chunk unique rows", "Chunk class", "Chunks total", "Chunks remaining"
      )
    )
    results <- lapply(chunk_info,process_chunk)
  }
  
  vebcat("File chunking protocol completed", color = "funSuccess")
}

#####################
# Chunk loadable dt
#####################

chunk_loaded_df <- function(dt, approach = "precautionary", chunk.name = "species", chunk.column, chunk.dir, iterations = NULL, verbose = FALSE) {
  vebcat("Initiating loaded dataframe chunking protocol.", color = "funInit")
  # File to store the last iteration number
  last_iteration_file <- paste0(chunk.dir, "/loaded-iteration.txt")
  
  n_total <- length(unique(dt[[chunk.column]]))
  
  # Check if the last iteration file exists
  if (file.exists(last_iteration_file)) {
    if (is.null(iterations)) {
      iterations <- as.integer(readLines(last_iteration_file))
      
      if (iterations >= n_total) {
        return(vebcat("Data already chunked", color = "funSuccess"))
      }
    }
  }
  
  dt <- select_species_approach(
    dt = dt,
    approach = approach,
    col.name = chunk.column,
    custom.list = config$species$taxonRank_infraEpithets,
    verbose = verbose
  )
  
  vebprint(unique(dt[[chunk.column]], veb = verbose, "Selected species:"))
  
  vebcat("Sorting data.", veb = verbose)
  data <- dt[order(dt[[chunk.column]]), ]
  
  vebcat("Splitting data into", highcat(chunk.column), "lists.", veb = verbose)
  data_list <- split(data, data[[chunk.column]])
  
  n_total <- length(data_list)
  
  create_dir_if(chunk.dir)
  create_dir_if(paste0(chunk.dir, "/", chunk.name))
  
  err_log <- paste0(chunk.dir, "/loaded-error.txt")
  create_file_if(err_log)
  
  # If iterations is not provided, run all iterations
  if (is.null(iterations)) {
    iterations <- seq_along(data_list)
  }
  
  catn("Chunking dataframe into files")
  cat(sprintf("%10s | %16s | %8s \n", paste0("n_", chunk.name), paste0("total n_", chunk.name), "Remaining"))
  for (i in iterations) {
    tryCatch(
      {
        spec <- names(data_list)[i]
        cat(sprintf("\r%10.0f | %16.0f | %8.0f", i, n_total, n_total - i))
        flush.console()
        
        # Replace spaces with config$species$file_separator in x
        spec <- gsub(" ", config$species$file_separator, spec)
        
        # Remove trailing config$species$file_separator, if any
        spec <- sub(paste0(config$species$file_separator, "$"), "", spec)
        
        # Create the filename
        file_name <- paste0(chunk.dir, "/", chunk.name, "/", spec, ".csv")
        
        if (!file.exists(file_name)) {
          fwrite(data_list[[i]], file_name, bom = T)
        }
        
        writeLines(as.character(i), last_iteration_file)
      },
      error = function(e) {
        try(err_log <- file(err_log, open = "at"))
        sink(err_log, append = T, type = "output")
        
        catn("Error in iteration", i, ":", e$message)
        
        sink(type = "output")
        close(err_log)
      }
    )
  }
  catn()
  
  vebcat("Loaded dataframe chunking protocol completed successfully.", color = "funSuccess")
}

################
# Clean chunks
################

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
