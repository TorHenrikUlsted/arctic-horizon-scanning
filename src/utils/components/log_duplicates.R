log_duplicates <- function(df, column, process, folder_name, file_name) {
  
  df <- as.data.frame(df)
  
  dup_rows <- which(duplicated(df[[column]]))
  dup_sp <- df[dup_rows, column]
  
  log_out_df <- data.frame(scientificName = dup_sp, rowNumber = dup_rows)
  
  dups_n <- nrow(log_out_df)
  
  log_path <- file.path("./outputs", process, "logs", folder_name)
  if (!dir.exists(log_path)) dir.create(log_path, recursive = T)
  
  output_log_path <- file.path(log_path, file_name)
  
  fwrite(log_out_df, file = output_log_path, bom = T)
  
  cat("Number of duplicates:", cc$aquamarine(dups_n), "\n")
  cat("More info found in: ", yellow(output_log_path), "\n")
}
