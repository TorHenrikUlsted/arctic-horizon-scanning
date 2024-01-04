source_all("./src/setup/components")

check_cpu_speed <- function(df.path, max.cores, sample.size = NULL, verbose, counter) {
  hostname <- system("hostname", intern = T)
  dir_path <- paste0("./outputs/setup/cpu/", hostname)
  
  if (file.exists(paste0(dir_path, "/estimated_time.txt"))) {
    estimated_time <- readLines(paste0(dir_path, "/estimated_time.txt"))
    
    estimated_time <- as.numeric(estimated_time)
  } else {
    df <- fread(df.path, sep = "\t")
    
    sample_size <- min(nrow(df), sample.size)
    
    # Sample a subset of the data
    subset <- df[sample(nrow(df), size = sample_size), ]
    
    start_time <- Sys.time()
    
    result <- check_syn_wfo(subset, column = colnames(subset), folder = dir_path, min(nrow(df), max.cores), verbose = verbose, counter = counter)
    
    end_time <- Sys.time()
    
    # Calculate the time it took to run your function on the subset
    time_taken <- difftime(end_time, start_time, units = "secs")
    
    # Estimate the time it would take to run the function on the full dataset
    estimated_time <- time_taken / sample_size
    
    writeLines(as.character(estimated_time), paste0(dir_path, "/estimated_time.txt"))
  }
  
  time.const <<- as.numeric(estimated_time)
  
  time.setup <<- as.numeric(30)
  
  cat("Time setup (sec):", cc$lightSteelBlue(time.setup), "\n")
  cat("time constant (sec):", cc$lightSteelBlue(estimated_time), "\n")
}

setup_raw_data <- function(column, test = NULL, max.cores, verbose, counter) {
  cat(blue("Setting up raw data. \n"))
  
  if (!is.null(test) && length(test) > 0) {
    test <- wrangle_test(test = test, column, verbose = verbose)
    
    checklist <- test
    
    } else {
    aba <- wrangle_aba(column, verbose = verbose)
    ambio <- wrangle_ambio(column, verbose = verbose)
    glonaf <- wrangle_glonaf(column, verbose = verbose)
    
    checklist <- c(aba, ambio, glonaf)
  }
  
  if (verbose) cat("dfs added to checklist: \n")
  if (verbose) print(names(checklist))
  
  checked_dfs <- syncheck_dfs(
    checklist, 
    column,
    out.dir = "./outputs/setup/wrangle", 
    max.cores = max.cores, 
    verbose = verbose, 
    counter = counter
    )
  
  if (verbose) cat("Combining no-matches. \n")
  
  if (all(sapply(checked_dfs, is.null))) {
    cat("All data frames already exist. \n")
    
  } else {
    combined_df <- data.frame()
    
    # Loop over each list in checked_dfs
    for(i in 1:length(checked_dfs)){
      if (!is.null(checked_dfs[[i]]$wfo_one_nomatch) && nrow(checked_dfs[[i]]$wfo_one_nomatch) > 0) {
        if (names(checked_dfs)[i] != "") {
          cat("Getting nomatches for:", cc$lightSteelBlue(names(checked_dfs)[i]), "\n")
          
          checked_dfs[[i]]$wfo_one_nomatch$dfOrigin <- names(checked_dfs)[i]
          cat("Adding an origin column with:", names(checked_dfs)[i], "\n")
        }
        
        # rbind the wfo_one_nomatch data frame from each list
        combined_df <- rbind(combined_df, checked_dfs[[i]]$wfo_one_nomatch)
      }
    }
    
    nomatch_combined <- combined_df[!duplicated(combined_df[[column]]), ]
    
    cat("Writing combined no-matches to:", yellow("./outputs/setup/wrangle"), "\n")
    
    if (nrow(nomatch_combined) > 0) {
      cat("There were", cc$lightSteelBlue(nrow(nomatch_combined)), "species without matches. \n")
      fwrite(nomatch_combined, "./outputs/setup/wrangle/combined-wfo-nomatch.csv", bom = T)
    } else {
      cat("There were", cc$lightSteelBlue(0), "species without any matches. \n")
    }
  }
  
  cat(cc$lightGreen("raw data setup completed successfully. \n"))
}

setup_sp <- function(test = T, big_test = F) {
  if (test == T) {
    sp_df <- run_test(big_test)
  } else {
    sp_df <- run_full(format)
    cat(blue("Initiating full run. \n"))
  }
  
  return(sp_df)
} # end setup
