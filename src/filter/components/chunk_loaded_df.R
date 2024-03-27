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