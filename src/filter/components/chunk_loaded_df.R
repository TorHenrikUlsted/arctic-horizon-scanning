chunk_loaded_df <- function(df, chunk.name = "output", chunk.column, chunk.dir, iterations = NULL, verbose = T) {
  cat(blue("Initiating loaded dataframe chunking protocol. \n"))
  
  if (is.vector(chunk.column) && length(chunk.column) > 1) {
    if (verbose) cat("Found vector more than 1 in length, combining columns:", cc$lightSteelBlue(chunk.column), ".\n")
    df$combined <- df[, do.call(paste, c(.SD, sep = " ")), .SDcols = chunk.column]
    chunk.column <- "combined"
  }
  
  if (verbose) cat("Sorting data. \n")
  data <- df[order(df[[chunk.column]]), ]
  
  if (verbose) cat("Splitting data into", chunk.column, "lists. \n")
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
      cat("This data has already been chunked. \n")
      return(cat(cc$lightGreen("Loaded dataframe chunking protocol completed successfully. \n")))
    }
  }
  
  # If iterations is not provided, run all iterations
  if(is.null(iterations)) {
    iterations <- seq_along(data_list)
  }
  
  cat("Chunking dataframe into files \n")
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
      
      cat("Error in iteration", i, ":", e$message, "\n")
      
      sink(type = "output")
      close(err_log)
    })
  }
  cat(cc$lightGreen("\nLoaded dataframe chunking protocol completed successfully. \n"))
}