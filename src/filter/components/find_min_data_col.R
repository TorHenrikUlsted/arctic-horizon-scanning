find_min_data_col <- function(file_path, verbose = T) {
  # Read the first 10 rows of the file
  df <- fread(file_path, nrows = 10)
  
  # Calculate the total character length of each column
  column_lengths <- sapply(df, function(x) sum(nchar(as.character(x))))
  
  # Find the column name with the least character length
  least_data_column <- names(column_lengths)[which.min(column_lengths)]
  
  if (verbose) cat("Column with the least data:", cc$lightSteelBlue(least_data_column), "\n")
  
  return(least_data_column)
}