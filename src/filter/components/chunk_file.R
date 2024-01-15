chunk_file <- function(file_path, chunk.name, chunk.column, chunk.dir, chunk.size = 1e6, iterations = NULL, verbose = T) {
  cat(blue("Initiating file chunking protocol. \n"))
  
  if (!file.exists(file_path)) {
    cat(red("Cannot find file_path, reading", file_path, "\n"))
  } else {
    cat("Using file", cc$lightSteelBlue(file_path), "\n")
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
    
    if (verbose) cat("The file has", cc$lightSteelBlue(total_rows), "total number of rows. \n")
    if (verbose) cat("chunk.size: ", cc$lightSteelBlue(chunk.size), "\n")
    total_chunks <- ceiling(total_rows / chunk.size)
    if (verbose) cat("Total chunks:", cc$lightSteelBlue(total_chunks), "\n")
    # chunks_per_core <- ceiling(total_chunks / cores.max)
    # if (verbose) cat("chunks per core:", cc$lightSteelBlue(chunks_per_core), "\n")
    
    iteration_file <- paste0(chunk.dir, "/file-iteration.txt")
    
    if (file.exists(iteration_file)) {
      i_start <- as.numeric(readLines(iteration_file))
    } else {
      create_file_if(iteration_file)
      i_start <- 1
    }
  }
  
  if (i_start > total_chunks) {
    cat(cc$lightCoral("The starting iteration number is above total_chunks. Stopping the chunking process.\n"))
    return()
  }
  
  if (is.vector(chunk.column) && length(chunk.column) > 1) if (verbose) cat("Combining columns:", cc$lightSteelBlue(chunk.column), ".\n")
  
  df_header <- fread(file_path, nrows = 0)
  print(df_header)
  
  # cl <- makeCluster(cores.max)
  # on.exit(closeAllConnections())
  # 
  # clusterEvalQ(cl, {
  #   library(data.table)
  # })
  
  cat("Writing out process to:", yellow(chunk_path), "\n")
  
  cat("Chunking file into files \n")
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
      cat("Error occurred when reading data file with fread. \n")
      cat(e$message, "\n")
    })

      new_column <- chunk.column
      tryCatch({
        
        if (is.vector(chunk.column) && length(chunk.column) > 1) {
          data$combined <- apply(data[, ..chunk.column, drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " "))
          data$combined <- trimws(data$combined)
          new_column <- "combined"
        }
      }, error = function(e) {
        cat("Error when combining columns. \n")
        cat(e$message)
        
        cat("Length unique combined column:", length(unique(data[[new_column]])), "\n")
        cat("Str new_column:\n")
        print(str(head(new_column, 2)))
      })
      
      # Order and split the data
      tryCatch({
        data <- data[order(data[[new_column]]), ]
        data_list <- split(data, data[[new_column]])
        # Remove list items with blank names
        data_list <- data_list[names(data_list) != ""]
      }, error = function(e) {
        cat("Error when ordering and splitting data_list. \n")
        cat(e$message)
        
        cat("Length of data_list: ", length(data_list), "\n")
        cat("Class: ", class(data_list), "\n")
        cat("List names:\n")
        print(head(names(data_list, 2)))
        cat("str:\n")
        print(str(head(data_list, 2)))
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
          #cat("length of x", x, "\n")
          file_name <- paste0(chunk.dir, "/", chunk.name, "/", gsub(" ", "-", x), ".csv")
          #cat("File name:", file_name, "\n")
          
          if(!file.exists(file_name)) {
            fwrite(data_list[[x]], file_name, bom = T)
          } else {
            fwrite(data_list[[x]], file_name, bom = T, append = T)
          }
        })
      }, error = function(e) {
        cat("Error when applying names and writing out files. \n")
        cat(e$message)
        cat("Class: ", class(data_list), "\n")
        cat("List names:\n")
        print(head(names(data_list, 2)))
        cat("str:\n")
        print(str(head(data_list, 2)))
      })
      
      rm(data)
      rm(data_list)
      
      # And then call the garbage collector
      invisible(gc())
      
    }, error = function(e) {
      try(err_log <- file(err_log, open = "at"))
      sink(err_log, append = T, type = "output")
      
      cat("Error in chunk", i, ":", e$message, "\n")

      sink(type = "output")
      close(err_log)
    })
    
    writeLines(as.character(i), iteration_file)
    
    i <- i + 1
  }
  
  cat(cc$lightGreen("\nFile chunking protocol completed successfully. \n"))
}