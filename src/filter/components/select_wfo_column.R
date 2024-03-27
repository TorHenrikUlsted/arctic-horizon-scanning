select_wfo_column <- function(filepath, col.unique, col.select = NULL, col.combine = NULL, pattern = "*.csv", verbose = F) {
  
  catn("Selecting WFO column")
  
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
      vebcat("\nFound vector more than 1 in length, combining columns:", highcat(col.unique), veb = verbose)
      df$refinedScientificName <- apply(df[, ..col.unique, drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " "))
      df$refinedScientificName <- trimws(df$refinedScientificName)
      col.unique <- "refinedScientificName"
    }
    
    if (!is.null(col.select)) {
      vebcat("Selecting the", highcat(col.select), "column(s) and using", highcat(col.unique), "as the unique column.", veb = verbose)
      df_sel <- df %>% 
        select(all_of(c(col.select, col.unique)))
    } else {
      vebcat("Using the", col.unique, "as unique column.", veb = verbose)
      df_sel <- df %>% 
        select(all_of(col.unique))
    }
    
    df_uniq <- df_sel[!duplicated(df_sel[[col.unique]]), ]
    
    n_orig <- nrow(df_sel)
    n_uniq <- nrow(df_uniq)
    catn("\nList:", highcat(sub("-wfo-one.csv$", "", basename(file))))
    cat(sprintf("%-10s | %s \n", "n_species", "unique(n_species)"))
    cat(highcat(sprintf("%-10d | %d \n", n_orig, n_uniq)))
    catn()
    
    return(df_uniq)
  })
  
  names(df_list) <- sub("-wfo-one.csv$", "", basename(csv_files))
  names(df_list) <- gsub("-", "_", names(df_list))
  
  
  return(df_list = df_list)
}