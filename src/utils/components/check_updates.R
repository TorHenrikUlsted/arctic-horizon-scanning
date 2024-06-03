get_dependencies <- function(package.name, suggested.install = TRUE, recursive = TRUE) {
  depends <- packageDescription(package.name)$Depends
  
  
  find_dependencies <- function(x) {
    if (is.null(x)) return()
    
    x <- gsub("\n", "", x)
    
    split_strings <- c()
    
    if (length(strsplit(x, ", ")[[1]]) > 0) {
      split_x <- strsplit(x, ", ")[[1]]
      
      for (sx in split_x) {
        split_strings <- c(split_strings, strsplit(sx, " ")[[1]][[1]])
      }
      
    } else {
      split_strings <- c(split_strings, strsplit(x, " ")[[1]][[1]])  
    }
    
    split_strings <- gsub("\\(>=", "", split_strings)
    
    split_strings <- split_strings[!split_strings %in% "R"]
    
    return(split_strings)
  }
  
  dependencies <- find_dependencies(depends)
  
  if (suggested.install) {
    suggested <- packageDescription(package.name)$Suggests
    
    suggested_deps <- find_dependencies(suggested)
    
  }
  
  d_recursive <- c()
  s_recursive <- c()
  
  if (recursive) {
    
    for (x in suggested_deps) {
      cat("Checking for dependencies of", x, "\n")
      
      depends <- packageDescription(x)$Depends
      
      fd <- find_dependencies(depends)
      
      if (!is.null(fd)) d_recursive <- c(d_recursive, fd)
      
      suggests <- packageDescription(x)$Suggests
      
      fs <- find_dependencies(suggests)
      
      
      if (!is.null(fs))  s_recursive <- c(s_recursive, fs)
      
      
    }
  }
  
  combined <- c(dependencies, d_recursive, suggested_deps, s_recursive)
  
  combined <- combined[!duplicated(combined)]
  
  return(combined)
}

check_updates <- function(pkgs) {
  message("Setting up packages.")
  
  restart_needed <- FALSE
  
  message("Checking for outdated packages.")
  
  # Check for outdated packages without considering dependencies
  outdated <- intersect(old.packages(
    lib.loc = .libPaths()[1],
    type = getOption("pkgType"))[, "Package"], 
    pkgs
  )
  
  # Inform about outdated source packages
  if (length(outdated) > 0) {
    cat("Outdated packages found: ", paste(outdated, collapse = ", "), "\n")
    
    # Ask the user if they want to update the outdated packages
    update <- readline(prompt = "Do you want to update the outdated packages? [y/n] ")
    
    if (tolower(update) == "y") {
      cat("Updating outdated packages...\n")
      
      # Try to install binary packages
      for (pkg in outdated) {
        cat("Updating", pkg, "\n")
        dependencies <- get_dependencies(pkg, recursive = FALSE)
        
        depend_loc <- find.package(dependencies, quiet = TRUE)
        
        for (i in depend_loc) {
          if (file.access(i, 2) == 0) {
            file_access <- FALSE
          }
        }
        
        # Check if the directory is writable
        if (file_access) {
          # If the directory is writable, install the package in the current library path
          update.packages(
            oldPkgs = pkg, 
            ask = FALSE, 
            dependencies = TRUE, 
            repos = options()$repos, 
            type = getOption("pkgType")
          )
        } else {
          # If the directory is not writable, print a message and install the package in a personal library
          cat("The library path for some dependencies in the package ", pkg, " is not writable. Installing the dependencies to the personal library.\n")
          
          install.packages(pkg)
          
          # for (i in dependencies) {
          #   print(i)
          # 
          #   suppressWarnings(
          #     install.packages(
          #       pkgs = i,
          #       lib = .libPaths()[1],
          #       ask = FALSE,
          #       dependencies = TRUE,
          #       repos = options()$repos,
          #       type = getOption("pkgType")
          #     )
          #   )
          # }
        }
      }
      
      restart_needed <- TRUE
    } else {
      cat("Update process cancelled by user.\n")
    }
  } else {
    cat("All specified packages are up to date.\n")
  }
  
  cat("Installing and loading packages. \n")
  
  # Install packages if necessary and load them
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
        # If the directory is writable, install the package in the current library path
        install.packages(
          pkg, 
          lib = .libPaths()[1], 
          dependencies = TRUE, 
          repos = options()$repos, 
          type = getOption("pkgType")
        )
      
      
      library(pkg, character.only = TRUE)
      cat("Package ", pkg, " loaded successfully.\n")
      restart_needed <- TRUE
    }
  }
  
  
  return(restart_needed)
}
