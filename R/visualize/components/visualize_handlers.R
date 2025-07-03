check_visualization_files <- function(files, dir, vis.title, vis.save.device, verbose = FALSE) {
  # Convert single filename to list for consistency
  if (is.character(files) && length(files) == 1) {
    files <- list(files)
  }
  
  # Add file extension if not present
  files <- sapply(files, function(x) {
    if (!grepl(paste0("\\.", vis.save.device, "$"), x)) {
      paste0(x, ".", vis.save.device)
    } else {
      x
    }
  })
  
  if (dir.exists(dir)) {
    existing_files <- basename(list.files(dir, recursive = TRUE))
    missing_files <- setdiff(files, existing_files)
    
    if (verbose) {
      catn("Checking in directory:", dir)
      catn("Found", length(existing_files), "files")
    }
    
    return(list(
      complete = length(missing_files) == 0,
      missing = missing_files,
      existing = intersect(files, existing_files)
    ))
  }
  
  # If directory doesn't exist, everything is missing
  return(list(
    complete = FALSE,
    missing = files,
    existing = character(0)
  ))
}