# Function to get current memory usage
get_mem_usage <- function(type = "free", format = "b") {
  # Call free command and get the output
  
  mem_info <- strsplit(system("free -m", intern = TRUE)[2], "\\s+")[[1]]
  
  if (type == "free") {
    # Get free memory
    mem_usage <- as.numeric(mem_info[4])
  } else if (type == "total") {
    # Get total memory
    mem_usage <- as.numeric(mem_info[2])
  } else if (type == "used") {
    # Get used memory
    mem_usage <- as.numeric(mem_info[3])
  } else {
    stop("Invalid type. Choose 'free', 'total', or 'used'.")
  }
  
  # Transform to specified format
  if (format == "b") {
    mem_usage <- mem_usage * 1024^2
  } else if (format == "kb") {
    mem_usage <- mem_usage * 1024
  } else if (format == "mb") {
    mem_usage <- mem_usage
  } else if (format == "gb") {
    mem_usage <- mem_usage / 1024
  } else if (format == "tb") {
    mem_usage <- mem_usage / 1024^2
  } else {
    stop("Invalid format. Choose 'b', 'kb', 'mb', 'gb', or 'tb'.")
  }
  
  return(mem_usage)
}

start_mem_tracking <- function(file.out, stop_file) {
  init_val <- get_mem_usage(format = "gb")
  # Create a control object
  control <- new.env()
  
  # Start the tracking in a separate process
  control$pid <- parallel::mcparallel({
    # Initialize a vector to store memory usage over time
    mem_usage <- c()
    
    while(!file.exists(stop_file)) {
      # Check memory usage every second
      Sys.sleep(1)
      
      new_mem_usage <- get_mem_usage(format = "gb")
      
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