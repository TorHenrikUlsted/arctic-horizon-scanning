select_wfo_column <- function(read.dir, column, pattern = "*.csv", verbose = F) {
  
  csv_files <- list.files(path = read.dir, pattern = pattern, full.names = TRUE)
  
  if (verbose) cat("Selecting the", column, "column. \n")
  
  df_list <- lapply(csv_files, function(file) {
    df <- fread(file)
    
    df_sel <- df %>% 
      select(all_of(column))
    
    df_uniq <- unique(df_sel)
    
    n_orig <- nrow(df_sel)
    n_uniq <- nrow(df_uniq)
    
    if (verbose) {
      cat(sprintf("%-10s | %s \n", "n_species", "unique(n_species)"))
      cat(cc$lightSteelBlue(sprintf("%-10d | %d \n", n_orig, n_uniq)))
    }
    
    return(df_uniq)
  })
  
  names(df_list) <- sub("-wfo-one.csv$", "", basename(csv_files))
  names(df_list) <- gsub("-", "_", names(df_list))
  
  
  return(df_list)
}