#Handle all inputs and split into parts by "."
config$run <- list()
config_handler <- function(...) {
  call_args <- match.call(expand.dots = FALSE)$`...`
  args <- list(...)
  names(args) <- sapply(call_args, deparse)
  
  set_nested <- function(lst, keys, value) {
    if (length(keys) == 1) {
      lst[[keys]] <- value
    } else {
      if (is.null(lst[[keys[1]]])) {
        lst[[keys[1]]] <- list()
      }
      lst[[keys[1]]] <- set_nested(lst[[keys[1]]], keys[-1], value)
    }
    return(lst)
  }
  
  for (original_name in names(args)) {
    name <- original_name
    if (grepl(".", name, fixed = TRUE)) {
      name <- gsub(".", "$", name, fixed = TRUE)
    }
    name_parts <- strsplit(name, "$", fixed = TRUE)[[1]]
    
    config$run <<- set_nested(config$run, name_parts, args[[original_name]])
  }
}
