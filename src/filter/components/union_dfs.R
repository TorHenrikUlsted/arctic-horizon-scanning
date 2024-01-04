union_dfs <- function(df1, df2, verbose = F) {
  
  df1_name <- deparse(substitute(df1))
  df2_name <- deparse(substitute(df2))
  
  if (verbose) cat("Merging", cc$lightSteelBlue(df1_name), "and", cc$lightSteelBlue(df2_name), "\n")
  
  ## combine present lists and remove duplicates
  merged_df = dplyr::union(df1, df2)
  
  cat(sprintf("%13s | %13s | %9s \n", df1_name, df2_name, "merged_df"))
  cat(cc$lightSteelBlue(sprintf("%13d | %13d | %9d \n", nrow(df1), nrow(df2), nrow(merged_df))))
  
  cat("Duplicated species removed:", cc$lightSteelBlue(nrow(df1) + nrow(df2) - nrow(merged_df)), "\n")
  
  ## Run NA and distinct check
  if (any(is.na(merged_df)) == T) {
    cat(red("Some merged_df species are NA \n"))  
    
    merged_df <- na.omit(merged_df)
  } else {
    cat(green("No merged_df species are NA \n"))
  }
  
  if (any(duplicated(merged_df)) == T) {
    cat(red("Some merged_df species are duplicated \n"))  
    
    merged_df <- unique(merged_df)
  } else {
    cat(green("All merged_df species are unique \n"))
  }
  
  return(merged_df)
}