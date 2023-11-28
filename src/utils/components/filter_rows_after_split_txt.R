filter_rows_after_split_text <- function(df, col1, col2, split_text) {
  df %>%
    rowwise() %>%
    filter(!grepl(split_text, !!sym(col1)) |
      !grepl(split_text, !!sym(col2)) |
      (length(strsplit(!!sym(col1), split_text)[[1]]) > 1 &&
        length(strsplit(!!sym(col2), split_text)[[1]]) > 1 &&
        strsplit(!!sym(col1), split_text)[[1]][2] ==
          strsplit(!!sym(col2), split_text)[[1]][2]))
  
  return(df)
}
