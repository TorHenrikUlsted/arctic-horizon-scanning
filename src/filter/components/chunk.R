#####################
# Chunk from file
#####################

chunk_file <- function(file_path, chunk.name, chunk.column, chunk.dir, chunk.size = 1e6, iterations = NULL, verbose = T) {
  vebcat("Initiating file chunking protocol.", color = "funInit")
  
  if (!file.exists(file_path)) {
    vebcat("Cannot find file_path, reading", file_path, color = "nonFatalError")
  } else {
    catn("Using file", highcat(file_path))
  }
  
  create_dir_if(chunk.dir)
  chunk_path <- paste0(chunk.dir, "/", chunk.name)
  create_dir_if(chunk_path)
  
  # core_dir <- paste0(chunk.dir, "/nodes")
  # create_dir_if(core_dir)
  
  err_log <- paste0(chunk.dir, "/file-error.txt")
  create_file_if(err_log)
  
  # If iterations is not provided, calculate the total number of chunks
  if(is.null(iterations)) {
    col_min_data <- find_min_data_col(file_path)
    
    if (!file.exists(paste0(chunk.dir, "/total-rows.txt"))) {
      # Check the operating system
      if (Sys.info()["sysname"] == "Windows") {
        total_rows <- as.numeric(system2("findstr", args = c("/R", "/N", "^", file_path), stdout = TRUE, stderr = NULL))
      } else {  # for Unix-based systems like Linux and macOS
        total_rows <- as.numeric(system(paste("awk 'END {print NR}' ", file_path), intern = TRUE))
      }
      
      writeLines(as.character(total_rows), paste0(chunk.dir, "/total-rows.txt"))
    } else {
      total_rows <- as.numeric(readLines(paste0(chunk.dir, "/total-rows.txt")))
    }
    
    vebcat("The file has", highcat(total_rows), "total number of rows.", veb = verbose)
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
    vebcat("The starting iteration number is above total_chunks. Stopping the chunking process.", color = "nonFatalError")
    return()
  }
  
  if (is.vector(chunk.column) && length(chunk.column) > 1) {
    vebcat("Combining columns:", highcat(chunk.column), veb = verbose)
  }
  
  df_header <- fread(file_path, nrows = 0)
  print(df_header)
  
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
      "%6s | %14s | %20s |  %13s | %14s | %18s\n", 
      "Chunk", "Chunk n_rows", "Chunk n_uniq_rows", "Chunk class", "Chunks total", "Chunks remaining"
    )
  )
  
  
  i <- i_start
  while(TRUE) {
    tryCatch({
      
      if (i > total_chunks) {
        # Write the current iteration number to the file
        writeLines(as.character(i), iteration_file)
        break
      }
      
      # Read a chunk of the data
      tryCatch({
      data <- fread(file_path, skip = (i - 1) * chunk.size + 1, nrows = chunk.size, col.names = names(df_header), verbose = F)
    }, error = function(e) {
      catn("Error occurred when reading data file with fread.")
      catn(e$message)
    })

      new_column <- chunk.column
      tryCatch({
        
        if (is.vector(chunk.column) && length(chunk.column) > 1) {
          data$combined <- apply(data[, ..chunk.column, drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " "))
          data$combined <- trimws(data$combined)
          new_column <- "combined"
        }
      }, error = function(e) {
        catn("Error when combining columns.")
        catn(e$message)
        
        catn("Length unique combined column:", length(unique(data[[new_column]])))
        catn("Str new_column:")
        print(str(head(new_column, 2)))
      })
      
      # Order and split the data
      tryCatch({
        data <- data[order(data[[new_column]]), ]
        data_list <- split(data, data[[new_column]])
        # Remove list items with blank names
        data_list <- data_list[names(data_list) != ""]
      }, error = function(e) {
        catn("Error when ordering and splitting data_list.")
        catn(e$message)
        
        catn("Length of data_list: ", length(data_list))
        catn("Class: ", class(data_list))
        vebprint(head(names(data_list, 2)), text = "List names:")
        print(str(head(data_list, 2)), text = "str:")
      })
      
      cat(
        sprintf(
          "\r%6.0f | %14.0f | %20.0f | %13s | %14.0f | %18.0f",
          i, nrow(data), length(unique(data[[new_column]])), class(data_list), total_chunks, ceiling((total_rows - (i * chunk.size)) / chunk.size)
        )
      )
      flush.console()
      
      tryCatch({
        lapply(names(data_list), function(x) {
          #catn("length of x", x)
          file_name <- paste0(chunk.dir, "/", chunk.name, "/", gsub(" ", "-", x), ".csv")
          #catn("File name:", file_name)
          
          if(!file.exists(file_name)) {
            fwrite(data_list[[x]], file_name, bom = T)
          } else {
            fwrite(data_list[[x]], file_name, bom = T, append = T)
          }
        })
      }, error = function(e) {
        catn("Error when applying names and writing out files.")
        catn(e$message)
        vebprint(class(data_list), text = "Class:")
        vebprint(head(names(data_list, 3)), text = "List names:")
        vebprint(str(head(data_list, 3)), text = "str:")
      })
      
      rm(data)
      rm(data_list)
      
      # And then call the garbage collector
      invisible(gc())
      
    }, error = function(e) {
      try(err_log <- file(err_log, open = "at"))
      sink(err_log, append = T, type = "output")
      
      catn("Error in chunk", i, ":", e$message)

      sink(type = "output")
      close(err_log)
    }); catn()
    
    writeLines(as.character(i), iteration_file)
    
    i <- i + 1
  }
  
  vebcat("File chunking protocol completed successfully", veb = verbose)
}

#####################
# Chunk loadable df
#####################

chunk_loaded_df <- function(df, chunk.name = "output", chunk.column, chunk.dir, iterations = NULL, verbose = T) {
  vebcat("Initiating loaded dataframe chunking protocol.", color = "funInit")
  
  if (is.vector(chunk.column) && length(chunk.column) > 1) {
    vebcat("Found vector more than 1 in length, combining columns:", highcat(chunk.column), veb = verbose)
    df$combined <- df[, do.call(paste, c(.SD, sep = " ")), .SDcols = chunk.column]
    chunk.column <- "combined"
  }
  
  vebcat("Sorting data.", veb = verbose)
  data <- df[order(df[[chunk.column]]), ]
  
  vebcat("Splitting data into", highcat(chunk.column), "lists.", veb = verbose)
  data_list <- split(data, data[[chunk.column]])
  
  n_total <- length(data_list)
  
  create_dir_if(chunk.dir)
  create_dir_if(paste0(chunk.dir, "/", chunk.name))
  
  err_log <- paste0(chunk.dir, "/loaded-error.txt")
  if (!file.exists(err_log)) {
    file.create(err_log)
  } else {
    file.remove(err_log)
    file.create(err_log)
  }
  
  # File to store the last iteration number
  last_iteration_file <- paste0(chunk.dir, "/loaded-iteration.txt")
  
  # Check if the last iteration file exists
  if(file.exists(last_iteration_file)) {
    # Read the last iteration number from the file
    last_iteration <- as.integer(readLines(last_iteration_file))
    
    # Check if the last iteration is the same as the total number of iterations
    if(last_iteration >= n_total) {
      catn("This data has already been chunked.")
      return(vebcat("Loaded dataframe chunking protocol completed successfully.", color = "funSuccess"))
    }
  }
  
  # If iterations is not provided, run all iterations
  if(is.null(iterations)) {
    iterations <- seq_along(data_list)
  }
  
  catn("Chunking dataframe into files")
  cat(sprintf("%8s | %16s | %8s \n", paste0("n_", chunk.name), paste0("total n_", chunk.name), "Remaining"))
  for(i in iterations) {
    tryCatch({
      x <- names(data_list)[i]
      cat(sprintf("\r%8.0f | %16.0f | %8.0f", i, n_total, n_total - i))
      flush.console()
      
      # Replace spaces with hyphens in x
      x <- gsub(" ", "-", x)
      
      # Remove trailing hyphen, if any
      x <- sub("-$", "", x)
      
      # Create the filename
      file_name <- paste0(chunk.dir, "/", chunk.name, "/", x, ".csv")
      
      if(!file.exists(file_name)) {
        fwrite(data_list[[i]], file_name, bom = T)
      }
      
      writeLines(as.character(i), last_iteration_file)
    }, error = function(e) {
      try(err_log <- file(err_log, open = "at"))
      sink(err_log, append = T, type = "output")
      
      catn("Error in iteration", i, ":", e$message)
      
      sink(type = "output")
      close(err_log)
    })
  }; catn()
  
  vebcat("Loaded dataframe chunking protocol completed successfully.", color = "funSuccess")
}

################
# Clean chunks
################

clean_chunks <- function(chunk.name, chunk.column, chunk.dir, sp_w_keys, iterations = NULL, verbose = F) {
  vebcat("Initiating chunk cleaning protocol.", color = "funInit")
  
  # Get the list of all chunk files
  chunk_files <- list.files(paste0(chunk.dir, "/", chunk.name), pattern = "*.csv", full.names = TRUE)
  
  if(is.null(iterations)) {
    iterations <- seq_along(chunk_files)
  }
  
  err_log <- paste0(chunk.dir, "/cleaner-error.txt")
  if (!file.exists(err_log)) file.create(err_log) else {
    file.remove(err_log)
    file.create(err_log)
  }
  
  # Total number of chunk files
  n_total <- length(iterations)
  
  # File to store the last iteration number
  last_iteration_file <- paste0(chunk.dir, "/cleaner-iteration.txt")
  
  # Check if the last iteration file exists
  if(file.exists(last_iteration_file)) {
    # Read the last iteration number from the file
    last_iteration <- as.integer(readLines(last_iteration_file))
    
    # Check if the last iteration is the same as the total number of iterations
    if(last_iteration >= n_total) {
      catn("This data has already been cleaned.")
      return(vebcat("Chunk cleaning protocol completed successfully.", color = "funSuccess"))
    }
  }
  
  catn("Cleaning chunk files")
  cat(sprintf("%5s | %8s | %16s | %8s \n", "File", "n_data", "Total n_data", "Remaining"))
  
  j <- 1
  # Loop over each chunk file
  for(i in iterations) {
    tryCatch({
      # Progress indicator
      cat(sprintf("\r%5.0f | %8.0f | %16.0f | %8.0f", j, i, n_total, n_total - i))
      
      # Read the chunk file
      chunk_data <- fread(chunk_files[i])
      
      # Filter rows where scientificName matches with the scientificName of sp_w_keys
      cleaned_data <- chunk_data[combined %in% sp_w_keys$refinedScientificName, ]
      
      # Check if cleaned_data is empty
      if (nrow(cleaned_data) == 0) {
        # Remove the file
        file.remove(chunk_files[i])
      } else {
        # Write the cleaned data back to the file
        fwrite(cleaned_data, chunk_files[i], bom = T)
      }
      
      writeLines(as.character(i), last_iteration_file)
      
      j + 1
    }, error = function(e) {
      try(err_log <- file(err_log, open = "at"))
      sink(err_log, type = "output")
      
      catn("Error in iteration", i, ":", e$message)
      
      sink(type = "output")
      close(err_log)
    })
  }; catn()
  
  vebcat("Chunk cleaning protocol completed successfully.", color = "funSuccess")
}
