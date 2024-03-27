lock <- function(lock.dir, lock.n = 1) {
  # Check if the path is a directory or a file
  if (!file.info(lock.dir)$isdir) stop("Input string not a directory.")
  
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
      
      return(full_path)
    }
  }
  
  # If no lock file could be created, return an error
  cat("No available lock files. Increase lock.n or remove existing lock files.")
}

  


is.locked <- function(lock.dir, lock.n = 1) {
  locked_n <- length(list.files(lock.dir))
  
  if (locked_n < lock.n) return(FALSE) else return(TRUE)
}


unlock <- function(lock.object) {
  tryCatch({
    file.remove(lock.object)
  }, error = function(e) {
    stop("Failed to unlock ", lock.object, " with error: ", e$message)
  })
}
