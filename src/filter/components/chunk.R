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

  create_dir_if(chunk.dir)
  chunk_path <- paste0(chunk.dir, "/", chunk.name)
  create_dir_if(chunk_path)

  # core_dir <- paste0(chunk.dir, "/nodes")
  # create_dir_if(core_dir)

  err_log <- paste0(chunk.dir, "/file-error.txt")
  create_file_if(err_log)

  # If iterations is not provided, calculate the total number of chunks
  if (is.null(iterations)) {
    col_min_data <- find_min_data_col(file.path)

    if (!file.exists(paste0(chunk.dir, "/total-rows.txt"))) {
      total_rows <- system_calc_rows(file.path)

      writeLines(as.character(total_rows), paste0(chunk.dir, "/total-rows.txt"))
    } else {
      total_rows <- as.numeric(readLines(paste0(chunk.dir, "/total-rows.txt")))
    }

    vebcat("The file has", highcat(total_rows), "total number of rows.")
    vebcat("chunk.size: ", highcat(chunk.size), veb = verbose)
    total_chunks <- ceiling(total_rows / chunk.size)
    vebcat("Total chunks:", highcat(total_chunks), veb = verbose)
    # chunks_per_core <- ceiling(total_chunks / cores.max)
    # vebcat("chunks per core:", highcat(chunks_per_core), veb = verbose)

    iteration_file <- paste0(chunk.dir, "/file-iteration.txt")

    if (file.exists(iteration_file)) {
      i_start <- as.numeric(readLines(iteration_file))
    } else {
      create_file_if(iteration_file)
      i_start <- 1
    }
  }

  if (i_start > total_chunks) {
    return(vebcat("Data already chunked", color = "funSuccess"))
  }

  if (is.vector(chunk.column) && length(chunk.column) > 1) {
    vebcat("Combining columns:", highcat(chunk.column), veb = verbose)
  }

  df_header <- fread(file.path, nrows = 0)
  vebprint(df_header, verbose, "data table header:")

  # cl <- makeCluster(cores.max)
  # on.exit(closeAllConnections())
  #
  # clusterEvalQ(cl, {
  #   library(data.table)
  # })

  catn("Writing out process to:", colcat(chunk_path, color = "output"))

  catn("Chunking file into files \n")
  cat(
    sprintf(
      "%7s | %14s | %19s |  %13s | %13s | %17s\n",
      "Chunk", "Chunk n_rows", "Chunk n_uniq_rows", "Chunk class", "Chunks total", "Chunks remaining"
    )
  )

  i <- i_start
  while (TRUE) {
    tryCatch(
      {
        if (i > total_chunks) {
          # Write the current iteration number to the file
          writeLines(as.character(i), iteration_file)
          break
        }

        # Read a chunk of the data
        tryCatch(
          {
            data <- fread(file.path, skip = (i - 1) * chunk.size + 1, nrows = chunk.size, col.names = names(df_header), verbose = F)
          },
          error = function(e) {
            catn("Error occurred when reading data file with fread.")
            catn(e$message)
          }
        )

        new_column <- chunk.column

        tryCatch(
          {
           data <- select_species_approach(
              dt = data,
              approach = approach,
              col.name = chunk.column,
              custom.list = config$species$taxonRank_infraEpithet,
              verbose = verbose
            )
          },
          error = function(e) {
            catn("Error when selecting spcies approach.")
            catn("Approach:", approach)
            catn("Colmn name input:", chunk.column)
            vebprint(data, text = "Input data:")
            catn("Error message:")
            catn(e$message)

            catn("Length unique cleaned species name:", length(unique(data[[chunk.column]])))
            catn("Str chunk.column:")
            print(str(head(chunk.column, 2)))
          }
        )
        
        # Order and split the data
        tryCatch(
          {
            data <- data[order(data[[chunk.column]]), ]
            data_list <- split(data, data[[chunk.column]])
            # Remove list items with blank names
            data_list <- data_list[names(data_list) != ""]
          },
          error = function(e) {
            catn("Error when ordering and splitting data_list.")
            catn(e$message)
            
            catn("Length of data_list: ", length(data_list))
            catn("Class: ", class(data_list))
            vebprint(head(names(data_list, 2)), text = "List names:")
            print(str(head(data_list, 2)), text = "str:")
          }
        )
        
        cat(
          sprintf(
            "\r%7.0f | %14.0f | %19.0f | %13s | %13.0f | %17.0f",
            i, nrow(data), length(unique(data[[chunk.column]])), class(data_list), total_chunks, ceiling((total_rows - (i * chunk.size)) / chunk.size)
          )
        )
        flush.console()
        
        tryCatch(
          {
            lapply(names(data_list), function(x) {
              # Replace spaces with config$species$file_separator in x
              x <- gsub(" ", config$species$file_separator, x)
              
              # Remove trailing config$species$file_separator, if any
              x <- sub(paste0(config$species$file_separator, "$"), "", x)
              
              # Create the filename
              file_name <- paste0(chunk.dir, "/", chunk.name, "/", x, ".csv")
              
              if (!file.exists(file_name)) {
                fwrite(data_list[[x]], file_name, bom = T)
              } else {
                fwrite(data_list[[x]], file_name, bom = T, append = T)
              }
            })
          },
          error = function(e) {
            catn("Error when applying names and writing out files.")
            catn(e$message)
            vebprint(class(data_list), text = "Class:")
            vebprint(head(names(data_list, 3)), text = "List names:")
            vebprint(str(head(data_list, 3)), text = "str:")
          }
        )
        
        rm(c(data, data_list))
        invisible(gc())
      },
      error = function(e) {
        try(err_log <- file(err_log, open = "at"))
        sink(err_log, append = T, type = "output")
        
        catn("Error in chunk", i, ":", e$message)
        
        sink(type = "output")
        close(err_log)
      }
    )
    catn()
    
    writeLines(as.character(i), iteration_file)
    
    i <- i + 1
  }
  
  vebcat("File chunking protocol completed successfully", veb = verbose)
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
