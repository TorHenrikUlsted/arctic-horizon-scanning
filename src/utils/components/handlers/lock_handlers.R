lock <- function(lock.dir, lock.n = 1, lock.name = NULL) {
  # Check if the path is a directory or a file
  if (!file.info(lock.dir)$isdir) return(vebcat("Input string not a directory.", color = "nonFatalError"))
  
  # Check which lock files exist
  lock_exists <- rep(FALSE, lock.n)
  for (i in 1:lock.n) {
    full_path <- paste0(lock.dir, "/lock-", i, ".txt")
    if (file.exists(full_path)) {
      lock_exists[i] <- TRUE
    }
  }
  
  # Find the first lock file that does not exist and create it
  for (i in 1:lock.n) {
    if (!lock_exists[i]) {
      full_path <- paste0(lock.dir, "/lock-", i, ".txt")
      file.create(full_path)
      if (!is.null(lock.name)) {
        try(lock_con <- file(full_path, "w"))
        writeLines(as.character(lock.name), lock_con)
        close(lock_con)
      }
      
      return(full_path)
    }
  }
  
  # If no lock file could be created, return an error
  catn("No available lock files. Increase lock.n or remove existing lock files.")
}

  


is.locked <- function(lock.dir, lock.n = 1) {
  locked_n <- length(list.files(lock.dir))
  
  if (locked_n < lock.n) return(FALSE) else return(TRUE)
}


unlock <- function(lock.object) {
  if (is.null(lock.object)) {
    catn("Could not find lock:", lock.object, "continuing without unlocking.")
    return(invisible())
  }
  
  tryCatch({
    if (file.exists(lock.object)) {
      file.remove(lock.object)
    } else {
      return(invisible())
    }
  }, error = function(e) {
    return(vebcat("Failed to unlock ", lock.object, " with error: ", e$message, color = "nonFatalError"))
  })
}