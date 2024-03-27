union_dfs <- function(df1, df2, verbose = F) {
  
  catn("Combining data tables using union.")
  
  df1_name <- deparse(substitute(df1))
  df2_name <- deparse(substitute(df2))
  
  vebcat("Merging", highcat(df1_name), "and", highcat(df2_name), veb = verbose)
  
  ## combine present lists and remove duplicates
  merged_df = dplyr::union(df1, df2)
  
  cat(sprintf("%13s | %13s | %9s \n", df1_name, df2_name, "merged_df"))
  cat(highcat(sprintf("%13d | %13d | %9d \n", nrow(df1), nrow(df2), nrow(merged_df))))
  
  catn("Duplicated species removed:", highcat(nrow(df1) + nrow(df2) - nrow(merged_df)))
  
  ## Run NA and distinct check
  if (any(is.na(merged_df)) == T) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    
    merged_df <- na.omit(merged_df)
  }
  
  if (any(duplicated(merged_df)) == T) {
    vebcat("Some merged_df species are duplicated.", color = "nonFatalError")
    
    merged_df <- unique(merged_df)
  }
  
  return(merged_df)
}