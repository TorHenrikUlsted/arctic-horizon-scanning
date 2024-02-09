# Set encoding to UTF8
set_df_utf8 <- function(df) {
  
  for (name in names(df)[sapply(df, is.character)]) {
    df[[name]] <- enc2utf8(df[[name]])
  }
  
  return(df)
}