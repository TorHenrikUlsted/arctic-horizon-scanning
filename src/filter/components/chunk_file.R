chunk_file <- function(file_path, chunk.column, chunk.out, chunk.size = 1e6, iterations = NULL, verbose = T) {
  cat(blue("Initiating file chunking protocol. \n"))
  
  if (is.vector(chunk.column) && length(chunk.column) > 1) {
    if (verbose) cat("Found vector more than 1 in length, combining columns:", cc$lightSteelBlue(chunk.column), ".\n")
   df$combined <- df[, do.call(paste, c(.SD, sep = " ")), .SDcols = chunk.column]
    chunk.column <- "combined"
  }
  
  if (!file.exists(file_path)) {
    cat(red("Cannot find file_path, reading", file_path, "\n"))
  } else {
    cat("Using file", cc$lightSteelBlue(file_path), "\n")
  }
  
  create_dir_if(chunk.out)
  create_dir_if(paste0(chunk.out, "/", chunk.column))
  
  cat("Writing out files to:", yellow(paste0(chunk.out, "/", chunk.column)), "\n")
  
  err_log <- paste0(chunk.out, "/file-error.txt")
  if (!file.exists(err_log)) {
    file.create(err_log)
  } else {
    file.remove(err_log)
    file.create(err_log)
  }
  
  # File to store the iteration number
  iteration_file <- paste0(chunk.out, "/last-iteration.txt")
  
  # If iterations is not provided, calculate the total number of chunks
  if(is.null(iterations)) {
    total_rows <- nrow(fread(file_path, header = T))
    iterations <- seq_len(ceiling(total_rows / chunk.size))
  }
  
  cat("\nChunking file into files \n")
  cat(sprintf("%8s | %16s | %8s \n", "n_data", "Total n_data", "Remaining"))
  
  for(i in iterations) {
    tryCatch({
      cat(sprintf("\r%8.0f | %16.0f | %8.0f", i, length(iterations), length(iterations) - i))
      flush.console()
      
      # Read a chunk of the data
      data <- fread(file_path, skip = (i - 1) * chunk.size + 1, nrows = chunk.size)
      
      # Order and split the data
      data <- data[order(data[[chunk.column]]), ]
      data_list <- split(data, data[[chunk.column]])
      

      lapply(names(data_list), function(x) {
        file_name <- paste0(chunk.out, "/", chunk.column, "/", gsub(" ", "-", x), ".csv")
        
        if(!file.exists(file_name)) {
          fwrite(data_list[[x]], file_name, bom = T)
        } else {
          fwrite(data_list[[x]], file_name, bom = T, append = T)
        }
      })
      
      # Write the current iteration number to the file
      writeLines(as.character(i), iteration_file)
    }, error = function(e) {
      try(err_log <- file(err_log, open = "at"))
      sink(err_log, append = T, type = "output")
      
      cat("Error in chunk", i, ":", e$message, "\n")
      
      sink(type = "output")
      close(err_log)
    })
  }
  
  cat(cc$lightGreen("\nFile chunking protocol completed successfully. \n"))
}
