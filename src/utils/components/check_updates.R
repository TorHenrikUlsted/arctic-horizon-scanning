check_updates <- function(pkgs) {
  message("Setting up packages.")
  
  # Specify the log file path
  if (!dir.exists(paste0("./outputs/utils/logs"))) dir.create(paste0("./outputs/utils/logs"), recursive = T)
  
  log_file <- file.path("outputs/utils/logs/install_log.txt")
  
  # Check for missing packages
  missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing_pkgs) > 0) {
    cat("Missing packages found: ", paste(missing_pkgs, collapse = ", "), "\n")
    cat("Installing missing packages...\n")
    
    # Try to install binary packages
    sink(log_file)
    install.packages(missing_pkgs, dependencies = TRUE, repos = options()$repos, type = "binary")
    sink()
    search_string <- "type 'binary' is not supported"
    log <- readLines(log_file)
    is_present <- any(grepl(search_string, log, fixed = TRUE))
    
    # If installing binary packages failed, switch to source
    if (is_present) {
      cat("Failed to install binary packages, switching to source.\n")
      install.packages(missing_pkgs, dependencies = TRUE, repos = options()$repos, type = "source")
    } else {
      cat("Finished attempting to install binary packages.\n")
    }
    
    restart_needed <- TRUE
    cat("Restart needed: ", restart_needed, "\n")
  }
  
  message("Checking for outdated packages.")
  
  # Check for outdated packages without considering dependencies
  outdated_binary <- intersect(old.packages(type = "binary")[, "Package"], pkgs)
  outdated_source <- intersect(old.packages(type = "source")[, "Package"], pkgs)
  
  # Inform about outdated source packages
  if (length(outdated_source) > 0) {
    cat("Outdated source packages found: ", paste(outdated_source, collapse = ", "), "\n")
  }
  
  # Update outdated binary packages
  restart_needed <- FALSE
  if (length(outdated_binary) > 0) {
    cat("Outdated binary packages found: ", paste(outdated_binary, collapse = ", "), "\n")
    
    # Ask the user if they want to update the outdated packages
    update <- readline(prompt = "Do you want to update the outdated packages? [y/n] ")
    
    if (tolower(update) == "y") {
      cat("Updating outdated binary packages...\n")
      
      # Try to install binary packages
      sink(log_file)
      install.packages(outdated_binary, ask = FALSE, dependencies = TRUE, repos = options()$repos, type = "binary")
      sink()
      search_string <- "type 'binary' is not supported"
      log <- readLines(log_file)
      is_present <- any(grepl(search_string, log, fixed = TRUE))
      
      # If installing binary packages failed, switch to source
      if (is_present) {
        cat("Failed to install binary packages, switching to source.\n")
        install.packages(outdated_source, ask = FALSE, dependencies = TRUE, repos = options()$repos, type = "source")
      } else {
        cat("Finished attempting to install binary packages.\n")
      }
      
      restart_needed <- TRUE
    } else {
      cat("Binary packages not updated.\n")
    }
  } else {
    cat("All specified packages are up to date.\n")
  }
  
  cat("Installing and loading packages. \n")
  
  # Install packages if necessary and load them
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE, repos = options()$repos, type = "binary")
      library(pkg, character.only = TRUE)
      cat("Package ", pkg, " loaded successfully.\n")
      restart_needed <- TRUE
    }
  }
  
  return(restart_needed)
}
