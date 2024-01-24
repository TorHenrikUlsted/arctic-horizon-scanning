analyze_region_hv <- function(biovars, out.dir, name, method, verbose) {
  cat(blue("Analyzing cavm hypervolume. \n"))
  
  directory <- paste0(out.dir, "/", tolower(name), "/")
  
  create_dir_if(directory)
  
  if (!file.exists(paste0(directory, "hypervolume_", method, ".rds"))) {
    if (verbose) cat("File not found, initating hypervolume sequence. \n")
    matrix <- terra::values(biovars, na.rm = T)
    
    if (verbose) if (any(is.na(matrix))) cat(red("Some biovars region values are NA. \n")) else cat(green("No biovars region values are NA. \n"))
    
    if (verbose) cat("Matrix sample: \n")
    if (verbose) print(head(matrix, 3))
    
    hv <- hypervolume(matrix, name = name, method = method, verbose = verbose)
    
    # Save the hypervolume to a file
    saveRDS(hv, paste0(directory, "hypervolume_", method, ".rds"))
  } else {
    if (verbose) cat(name, "Hypervolume found, loading file. \n")
    # Load the hypervolume from the file
    hv <- readRDS(paste0(directory, "hypervolume_", method, ".rds"))
  }
  
  cat(cc$lightGreen(name, paste0("Hypervolume_", method), "Analysis Complete. \n"))
  return(hv)
}