format_mem_unit <- function(value, unit.in = "b", unit.out = "b") {
  
  units <- c("b", "kb", "mb", "gb", "tb")
  factors <- c(1, 1024, 1024^2, 1024^3, 1024^4)
  names(factors) <- units
  
  if (!(unit.in %in% units) | !(unit.out %in% units)) {
    return(vebcat("Invalid unit. Choose 'b', 'kb', 'mb', 'gb', or 'tb'.", color = "nonFatalError"))
  }
  
  value_bytes <- value * factors[unit.in]
  
  value_formatted <- value_bytes / factors[unit.out]
  
  names(value_formatted) <- unit.out
  
  return(value_formatted)
}

get_mem_usage <- function(type = "free", format = "b") {
  tryCatch({
    # Call free command and get the output
    mem_info <- system("free -m", intern = TRUE)
    
    # Check if we got the expected output
    if (length(mem_info) < 2) {
      stop("Unexpected output from 'free' command")
    }
    
    # Parse the memory information
    mem_values <- as.numeric(strsplit(mem_info[2], "\\s+")[[1]][-1])
    
    # Check if we parsed the expected number of values
    if (length(mem_values) < 3) {
      stop("Failed to parse memory information")
    }
    
    # Get the requested memory type
    mem_usage <- switch(type,
                        "free" = mem_values[3],
                        "total" = mem_values[1],
                        "used" = mem_values[2],
                        stop("Invalid type. Choose 'free', 'total', or 'used'."))
    
    # Format the memory value
    mem_usage <- format_mem_unit(mem_usage, unit.in = "mb", unit.out = format)
    
    return(unname(mem_usage))
  }, error = function(e) {
    warning(paste("Error in get_mem_usage:", e$message))
    return(NA)
  })
}

get_process_mem_use <- function(unit = "gb") {
  pid <- Sys.getpid()
  
  ram <- system2("ps", args=c("-p", pid, "-o", "rss="), stdout = TRUE)
  
  ram <- as.numeric(ram)
  
  ram <- format_mem_unit(ram, unit.in = "kb", unit.out = unit)
  
  ram <- round(ram, 3)
  
  catn("Current process RAM in", unit, highcat(unname(ram)))
  
  return(unname(ram))
}

log_message <- function(file.out, text, object = NULL) {
  try(con <- file(file.out, open = "a"))
  sink(con)
  if (!is.null(object)) {
    cat(paste0(Sys.time(), ":\n", text, "\n"))
    print(object)
  } else {
    cat(paste0(Sys.time(), ": ", text, "\n"))
  }
  sink()
  close(con)
}

start_mem_tracking <- function(file.out, file.stop) {
tryCatch({
  log_file <- file.path(dirname(file.out), "mem_tracking_log.txt")
  
  log_message(log_file, "Starting memory tracking")
  
  if(file.exists(file.stop)) {
    log_message("Removing existing stop file")
    file.remove(file.stop)
  }
  
  init_val <- get_mem_usage(type = "used", format = "gb")
  
  # Start the tracking in a separate R session
  system(paste0("Rscript -e \"
    source('./src/utils/components/memory_handlers.R')
    repeat {
      if (file.exists('", file.stop, "')) break
      new_mem_usage <- get_mem_usage(type = 'used', format = 'gb')
      if(file.exists('", file.out, "')) {
        existing_mem_usage <- as.numeric(readLines('", file.out, "'))
        max_mem_usage <- max(existing_mem_usage, new_mem_usage)
      } else {
        max_mem_usage <- new_mem_usage
      }
      writeLines(as.character(max_mem_usage), '", file.out, "')
      Sys.sleep(1)
    }
    \""), wait = FALSE)
  
}, error = function(e) {
  vebcat("Error when trying to start memory tracking", color = "fatalError")
  file.create(file.stop)
  stop(e)
})
  
  return(list(
    init.val = init_val,
    log = log_file,
    file.out = file.out,
    file.stop = file.stop
  ))
}

stop_mem_tracking <- function(control) {
  tryCatch({
    log_message(control$log, "Stopping memory tracking")
    # Create the stop file
    file.create(control$file.stop)
    
    # Give some time for the tracking process to finish
    Sys.sleep(2)
    
    # Read the final memory usage from the file
    max_mem_usage <- as.numeric(readLines(control$file.out))
    
    # Clean up files
    file.remove(control$file.out)
    if(file.exists(control$file.stop)) file.remove(control$file.stop)
    
    # Calculate peak memory usage
    peak_mem_usage <- max_mem_usage - control$init.val
    
    cat("Peak Memory Usage:", peak_mem_usage, "GB\n")
    
    # Write the peak memory usage to the output file
    writeLines(as.character(peak_mem_usage), control$file.out)
  }, error = function(e) {
    vebcat("Error when trying to stop memory tracking", color = "fatalError")
    stop(e)
  })
  
  return(peak_mem_usage)
}

# start_mem_tracking <- function(file.out, file.stop) {
#   if(file.exists(file.stop)) file.remove(file.stop)
#   create_file_if(file.out)
#   
#   init_val <- get_mem_usage(type = "used", format = "gb")
#   # Create a control object
#   control <- new.env()
#   
#   # Start the tracking in a separate process
#   control$pid <- parallel::mcparallel({
#     # Initialize a vector to store memory usage over time
#     mem_usage <- c()
#     
#     start_time <- Sys.time()
#     max_runtime <- 3600  # 1 hour, adjust as needed
#     
#     tryCatch({
#       while(!file.exists(file.stop) && (difftime(Sys.time(), start_time, units="secs") < max_runtime)) {
#         # Check memory usage every second
#         Sys.sleep(1)
#         
#         new_mem_usage <- get_mem_usage(type = "used", format = "gb")
#         
#         existing_mem_usage <- as.numeric(readLines(file.out))
#         
#         max_mem_usage <- max(c(existing_mem_usage, new_mem_usage), na.rm = TRUE)
#         
#         writeLines(as.character(max_mem_usage), file.out)
#       }
#     }, error = function(e) {
#       writeLines(paste("Error in child process:", e$message), file.stop)
#     }, finally = {
#       # Ensure the stop file is created even if there's an error
#       if(!file.exists(file.stop)) file.create(file.stop)
#     })
#   })
#   
#   return(list(
#     control = control, 
#     init.val = init_val,
#     file.out = file.out,
#     file.stop = file.stop
#   ))
# }
# 
# # Function to stop tracking memory usage
# stop_mem_tracking <- function(control) {
#   # Create the stop file
#   file.create(control$file.stop)
#   
#   # Wait for a short time to allow the child process to finish
#   Sys.sleep(2)
#   
#   # Try to collect the result, but don't wait indefinitely
#   result <- tryCatch({
#     parallel::mccollect(list(control$control$pid), wait = FALSE, timeout = 5)
#   }, error = function(e) {
#     message("Warning: Could not collect child process result. It may have already terminated.")
#     NULL
#   })
#   
#   # Read the memory usage over time from the file
#   mem_usage <- as.numeric(readLines(control$file.out))
#   
#   file.remove(control$file.out)
#   if(file.exists(control$file.stop)) file.remove(control$file.stop)
#   
#   # Get the peak memory usage
#   max_mem <- max(mem_usage)
#   
#   peak_mem_usage <- max_mem - control$init.val
#   
#   cat("Peak Memory Usage:\n")
#   print(peak_mem_usage)
#   
#   # Write the peak memory usage to an output file
#   writeLines(as.character(peak_mem_usage), control$file.out)
#   
#   return(peak_mem_usage)
# }

calc_num_cores <- function(ram.high, ram.low = 0, cores.total = detectCores(), verbose = FALSE) {
  # Ratio of high ram cores
  if (ram.low == 0) {
    high_core_ratio <- 1
    low_core_ratio <- 0
  } else {
    high_core_ratio <- (3/4)
    low_core_ratio <- (1/4)
  }
  
  # Calc ram the process will use
  high_load_cores <- floor(cores.total * high_core_ratio)
  low_load_cores <- floor(cores.total - high_load_cores)
  
  vebprint(cores.total, veb = verbose, "Total Cores:")
  vebprint(high_load_cores, veb = verbose, "High Load Cores:")
  vebprint(low_load_cores, veb = verbose, "Low Load Cores:")
  
  mem_limit_gb <- get_mem_usage("total", format = "gb")
  
  # Get the memory limit for high and low loads
  max_high_mem <- floor(mem_limit_gb * high_core_ratio)
  max_low_mem <- floor(mem_limit_gb * low_core_ratio)
  
  vebprint(max_high_mem, veb = verbose, "Max mem high:")
  vebprint(max_low_mem, veb = verbose, "Max mem low:")
  
  # Get the minimum of high load or the floor of max mem / peak ram
  max_cores_high <- min(high_load_cores, floor(max_high_mem / ram.high))
  if (low_load_cores == 0) {
    max_cores_low <- 0
  } else {
    max_cores_low <- min(low_load_cores, floor(max_low_mem / ram.low))
  }
  
  vebprint(max_cores_high, veb = verbose, "Max cores high:")
  vebprint(max_cores_low, veb = verbose, "Max cores low:")
  
  max_cores <- max_cores_high + max_cores_low
  
  max_cores <- min(max_cores, cores.total)
  
  vebprint(max_cores, veb = verbose, "max cores:")
  
  return(list(
    total = max_cores,
    high = max_cores_high,
    low = max_cores_low
  ))
}

get_disk_space <- function(directory = "/", units = "B") {
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    disk_space <- system2("cmd", args = "/c fsutil volume diskfree C:", stdout = TRUE)
    disk_space <- as.numeric(strsplit(disk_space[3], ": ")[[1]][2])
  } else if (os %in% c("Linux", "Darwin")) {
    disk_space <- system(paste0("df --output=avail -h ", directory, " | tail -n 1"), intern = TRUE)
    disk_space <- trimws(disk_space)  # Remove whitespace
    unit <- substr(disk_space, nchar(disk_space), nchar(disk_space))
    disk_space <- as.numeric(gsub(unit, "", disk_space))
    if (unit == "T") {
      disk_space <- disk_space * 1024  # Convert from terabytes to gigabytes
    }
  } else {
    stop("Unsupported operating system.")
  }
  
  conversion_factors <- c(B = 1, KB = 1024, MB = 1024^2, GB = 1024^3, TB = 1024^4)
  disk_space <- disk_space * 1e9 / conversion_factors[units]
  
  return(disk_space)
}




