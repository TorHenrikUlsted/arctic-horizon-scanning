check_system_speed <- function(df.path, test.name = "cpu-speed", sample.size = NULL, cores.max, fun, verbose = FALSE) {
  hostname <- system("hostname", intern = T)
  dir_path <- paste0("./outputs/setup/system/", hostname)
  filepath <- paste0(dir_path, "/", test.name, "-etc.txt")
  
  create_dir_if(dir_path)
  
  if (file.exists(filepath)) {
    estimated_time <- readLines(filepath)
    
    estimated_time <- as.numeric(estimated_time)
  } else {
    df <- fread(df.path, sep = "\t")
    
    sample_size <- min(nrow(df), sample.size)
    
    # Sample a subset of the data
    subset <- df[sample(nrow(df), size = sample_size), ]
    
    start_time <- Sys.time()
    
    result <- fun(
      subset = subset, 
      column = colnames(subset), 
      out.dir = dir_path, 
      cores.max = min(nrow(df), cores.max), 
      verbose = verbose
    )
    
    end_time <- Sys.time()
    
    # Calculate the time it took to run your function on the subset
    time_taken <- difftime(end_time, start_time, units = "secs")
    
    # Estimate the time it would take to run the function on the full dataset
    estimated_time <- time_taken / sample_size
    
    create_file_if(filepath)
    
    writeLines(as.character(estimated_time), filepath)
  }
  
  time_const <- as.numeric(estimated_time)
  
  catn("time constant (sec):", highcat(time_const))
  
  return(time_const)
}

wfo_speed <- function(subset, column, out.dir, cores.max, verbose = FALSE) {
  result <- check_syn_wfo(
    checklist = subset, 
    column = colnames(subset), 
    out.dir = out.dir, 
    cores.max = cores.max, 
    verbose = verbose, 
    counter = 1
  )
  
  return(result)
}
