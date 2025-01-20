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

get_disk_space <- function(directory = "/", units = "B") {
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    # For Windows, convert Unix-style path to Windows drive letter
    if (directory == "/") {
      directory <- "C:"  # Default to C drive
    }
    # Use PowerShell to get disk space in bytes
    disk_space <- system2("powershell", 
                          args = c("-command", 
                                   sprintf("(Get-PSDrive %s).Free",
                                           substr(directory, 1, 1))), 
                          stdout = TRUE)
    disk_space <- as.numeric(disk_space)
  } else if (os %in% c("Linux", "Darwin")) {
    disk_space <- system(paste0("df --output=avail -B1 ", directory, " | tail -n 1"), intern = TRUE)
    disk_space <- as.numeric(trimws(disk_space))
  } else {
    stop("Unsupported operating system.")
  }
  
  # Convert to requested units
  conversion_factors <- c(B = 1, KB = 1024, MB = 1024^2, GB = 1024^3, TB = 1024^4)
  disk_space <- disk_space / conversion_factors[units]
  
  return(disk_space)
}

get_mem_usage <- function(type = "free", format = "b") {
  os <- Sys.info()["sysname"]
  
  tryCatch({
    if (os == "Windows") {
      # Windows method using wmic
      mem_info <- system('wmic OS get FreePhysicalMemory,TotalVisibleMemorySize /Value', intern = TRUE)
      mem_info <- mem_info[mem_info != ""]
      total <- as.numeric(sub("TotalVisibleMemorySize=", "", mem_info[grep("TotalVisibleMemorySize", mem_info)]))
      free <- as.numeric(sub("FreePhysicalMemory=", "", mem_info[grep("FreePhysicalMemory", mem_info)]))
      used <- total - free
      
      # Convert from KB to requested format
      mem_values <- c(free = free, total = total, used = used)
      
    } else if (os %in% c("Linux", "Darwin")) {
      # Unix-like systems using free command
      mem_info <- system("free -b", intern = TRUE)
      if (length(mem_info) < 2) {
        # Try vmstat if free isn't available (e.g., on macOS)
        mem_info <- system("vm_stat", intern = TRUE)
        # Parse vm_stat output for macOS
        page_size <- 4096  # Default page size on macOS
        free <- as.numeric(gsub("[^0-9]", "", mem_info[grep("Pages free", mem_info)])) * page_size
        total <- as.numeric(gsub("[^0-9]", "", mem_info[grep("Pages active", mem_info)])) * page_size
        used <- total - free
      } else {
        # Parse free command output
        mem_values <- as.numeric(strsplit(mem_info[2], "\\s+")[[1]][-1])
        free <- mem_values[3]
        total <- mem_values[1]
        used <- mem_values[2]
      }
      mem_values <- c(free = free, total = total, used = used)
    } else {
      stop("Unsupported operating system")
    }
    
    # Get the requested memory type
    mem_usage <- switch(type,
                        "free" = mem_values["free"],
                        "total" = mem_values["total"],
                        "used" = mem_values["used"],
                        stop("Invalid type. Choose 'free', 'total', or 'used'."))
    
    # Format the memory value
    mem_usage <- format_mem_unit(mem_usage, unit.in = "b", unit.out = format)
    
    return(unname(mem_usage))
  }, error = function(e) {
    warning("Error getting memory usage: ", e$message, 
            "\nReturning NA - some features may be limited.")
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

track_memory <- function(fun, tracking = config$memory$tracking, identifier = NULL) {
  if (!tracking) return(fun)
  
  function(...) {
    fun_name <- if(!is.null(identifier)) identifier else deparse(substitute(fun))
    mem_start <- get_mem_usage("used", "gb")
    time_start <- Sys.time()
    
    # Run the function with cleanup
    result <- tryCatch({
      fun(...)
    }, finally = {
      mem_end <- get_mem_usage("used", "gb")
      time_end <- Sys.time()
      mem_diff <- mem_end - mem_start
      time_diff <- difftime(time_end, time_start, units = "mins")
      
      cat(sprintf("%s Memory: %0.2f GB, Time: %0.2f minutes\n", 
                  fun_name, mem_diff, as.numeric(time_diff)))
    })
    
    return(result)
  }
}

mem_check <- function(identifier = NULL, ram.use = NULL, custom_msg = NULL, interval = 60, verbose = FALSE) {
  tryCatch({
    mem_used_gb <- round(get_mem_usage(type = "used", format = "gb"), 2)
    mem_limit_gb <- round(config$memory$mem_limit / 1024^3, 2)
    
    # Static variable to track if message has been written for this identifier
    msg_written <- if (exists(".msg_written", envir = .GlobalEnv)) 
      get(".msg_written", envir = .GlobalEnv) 
    else 
      list()
    
    if (mem_used_gb >= mem_limit_gb) {
      if (!is.null(ram.use) && (!identifier %in% names(msg_written) || !msg_written[[identifier]])) {
        tryCatch({
          ram_con <- file(ram.use, open = "a")
          msg <- paste0(
            "RAM usage ", mem_used_gb, 
            " GB is above the maximum ", mem_limit_gb, " GB"
          )
          if (!is.null(custom_msg)) {
            msg <- paste0(msg, ". ", custom_msg)
          } else if (!is.null(identifier)) {
            msg <- paste0(msg, ". Waiting with ", identifier)
          }
          writeLines(msg, ram_con)
          close(ram_con)
          
          # Mark message as written for this identifier
          msg_written[[identifier]] <- TRUE
          assign(".msg_written", msg_written, envir = .GlobalEnv)
        }, error = function(e) {
          warning(paste("Failed to write to RAM usage file:", e$message))
        }, finally = {
          if(exists("ram_con")) try(close(ram_con))
        })
      }
      
      vebcat("Memory usage:", mem_used_gb, "GB exceeds limit:", mem_limit_gb, "GB", veb = verbose)
      invisible(gc(full = TRUE))
      
      if (interval > 0) {
        Sys.sleep(interval)
        Sys.sleep(runif(1, 0, 1))
      }
      
      return(TRUE)
    }
    
    # Reset message written status when memory is OK again
    if (!is.null(identifier) && identifier %in% names(msg_written)) {
      msg_written[[identifier]] <- FALSE
      assign(".msg_written", msg_written, envir = .GlobalEnv)
    }
    
    return(FALSE)
  }, error = function(e) {
    warning(paste("Error in memory check:", e$message))
    return(FALSE) # Return false on error to avoid infinite loops
  })
}

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
