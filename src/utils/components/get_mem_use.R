# Function to get current memory usage
get_mem_usage <- function() {
  # Call free command and get the output
  mem_usage <- as.numeric(strsplit(system("free -m", intern = TRUE)[2], "\\s+")[[1]][4])
  
  # Transform to bytes
  mem_usage <- mem_usage / 1024
  
  return(mem_usage)
}

start_mem_tracking <- function(file.out, stop_file) {
  init_val <- get_mem_usage()
  # Create a control object
  control <- new.env()
  
  # Start the tracking in a separate process
  control$pid <- parallel::mcparallel({
    # Initialize a vector to store memory usage over time
    mem_usage <- c()
    
    while(!file.exists(stop_file)) {
      # Check memory usage every second
      Sys.sleep(1)
      
      new_mem_usage <- get_mem_usage()
      
      existing_mem_usage <- as.numeric(readLines(file.out))
      
      min_mem_usage <- min(c(existing_mem_usage, new_mem_usage), na.rm = TRUE)
      
      writeLines(as.character(min_mem_usage), file.out)
    }
  })
  
  return(list(
    control = control, 
    init_val = init_val,
    file_out = file.out
  ))
}

# Function to stop tracking memory usage
stop_mem_tracking <- function(control, file.out, stop_file) {
  # Create the stop file
  file.create(stop_file)
  
  result <- parallel::mccollect(list(control$control$pid), wait = TRUE)
  
  # Read the memory usage over time from the file
  mem_usage <- as.numeric(readLines(control$file_out))
  
  # Get the peak memory usage
  min_mem <- min(mem_usage)
  
  print(min_mem)
  
  peak_mem_usage <- control$init_val - min_mem
  
  # Write the peak memory usage to an output file
  writeLines(as.character(peak_mem_usage), paste0(file.out))
  
  file.remove(control$file_out)
  file.remove(stop_file)
}