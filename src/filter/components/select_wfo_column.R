select_wfo_column <- function(filepath, col.unique, col.select = NULL, col.combine = NULL, pattern = "*.csv", verbose = F) {
  
  # Check if filepath is a directory or a specific file
  if (file.info(filepath)$isdir) {
    csv_files <- list.files(path = filepath, pattern = "*.csv", full.names = TRUE)
  } else {
    csv_files <- filepath
  }
  
  # make a list of data frames based on the different CSV files and also check for any "no matches" then add those to their own data frame.
  df_list <- lapply(csv_files, function(file) {
    df <- fread(file)
    
    if (is.vector(col.unique) && length(col.unique) > 1) {
      if (verbose) cat("\nFound vector more than 1 in length, combining columns:", cc$lightSteelBlue(col.unique), ".\n")
      df$refinedScientificName <- apply(df[, ..col.unique, drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " "))
      df$refinedScientificName <- trimws(df$refinedScientificName)
      col.unique <- "refinedScientificName"
    }
    
    if (!is.null(col.select)) {
      if (verbose) cat("Selecting the", cc$lightSteelBlue(col.select), "column(s) and using", cc$lightSteelBlue(col.unique), "as the unique column. \n")
      df_sel <- df %>% 
        select(all_of(c(col.select, col.unique)))
    } else {
      if (verbose) cat("Using the", col.unique, "as unique column. \n")
      df_sel <- df %>% 
        select(all_of(col.unique))
    }
    
    df_uniq <- df_sel[!duplicated(df_sel[[col.unique]]), ]
    
    n_orig <- nrow(df_sel)
    n_uniq <- nrow(df_uniq)
    cat("\nList:", cc$lightSteelBlue(sub("-wfo-one.csv$", "", basename(file))), "\n")
    cat(sprintf("%-10s | %s \n", "n_species", "unique(n_species)"))
    cat(cc$lightSteelBlue(sprintf("%-10d | %d \n", n_orig, n_uniq)))
    cat("\n")
    
    return(df_uniq)
  })
  
  names(df_list) <- sub("-wfo-one.csv$", "", basename(csv_files))
  names(df_list) <- gsub("-", "_", names(df_list))
  
  
  return(df_list = df_list)
}