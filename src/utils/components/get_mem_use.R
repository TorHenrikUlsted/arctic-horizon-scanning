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

start_mem_tracking <- function(file.out, file.stop) {
  if(file.exists(file.stop)) file.remove(file.stop)
  
  init_val <- get_mem_usage(type = "used", format = "gb")
  # Create a control object
  control <- new.env()
  
  # Start the tracking in a separate process
  control$pid <- parallel::mcparallel({
    # Initialize a vector to store memory usage over time
    mem_usage <- c()
    
    while(!file.exists(file.stop)) {
      # Check memory usage every second
      Sys.sleep(1)
      
      new_mem_usage <- get_mem_usage(type = "used", format = "gb")
      
      existing_mem_usage <- as.numeric(readLines(file.out))
      
      max_mem_usage <- max(c(existing_mem_usage, new_mem_usage), na.rm = TRUE)
      
      writeLines(as.character(max_mem_usage), file.out)
    }
  })
  
  return(list(
    control = control, 
    init.val = init_val,
    file.out = file.out
  ))
}

# Function to stop tracking memory usage
stop_mem_tracking <- function(control, file.stop) {
  # Create the stop file
  file.create(file.stop)
  
  result <- parallel::mccollect(list(control$control$pid), wait = TRUE)
  
  # Read the memory usage over time from the file
  mem_usage <- as.numeric(readLines(control$file.out))
  
  file.remove(control$file.out)
  
  # Get the peak memory usage
  max_mem <- max(mem_usage)
  
  peak_mem_usage <- max_mem - control$init.val
  
  catn("Peak Memory Usage:")
  print(peak_mem_usage)
  
  # Write the peak memory usage to an output file
  writeLines(as.character(peak_mem_usage), control$file.out)
  
  return(peak_mem_usage)
}
