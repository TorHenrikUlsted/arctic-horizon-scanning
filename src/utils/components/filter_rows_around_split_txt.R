filter_rows_around_split_text <- function(df, col1, col2, split_text) {
  df %>% filter(!grepl(split_text, !!sym(col1)) & !grepl(split_text, !!sym(col2)))

  return(df)
}
