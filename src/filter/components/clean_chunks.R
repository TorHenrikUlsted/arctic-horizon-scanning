clean_chunks <- function(chunk.name, chunk.column, chunk.dir, sp_w_keys, iterations = NULL, verbose = F) {
  cat(blue("Initiating chunk cleaning protocol. \n"))
  
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
      cat("This data has already been cleaned.\n")
      return(cat(cc$lightGreen("Chunk cleaning protocol completed successfully. \n")))
    }
  }
  
  cat("Cleaning chunk files \n")
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
      
      cat("Error in iteration", i, ":", e$message, "\n")
      
      sink(type = "output")
      close(err_log)
    })
  }
  
  cat(cc$lightGreen("\nChunk cleaning protocol completed successfully. \n"))
}
