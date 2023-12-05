source_all <- function(dir) {
  # Get a list of all .R files in the directory and its subdirectories
  r_files <- list.files(path = dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  
  # Source each file
  lapply(r_files, source)
  
  cat(length(r_files), "scripts sourced from", dir, "and its subdirectories.\n")
}
