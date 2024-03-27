print_function_args <- function() {
  # Get the function call
  fun_call <- sys.call(-1)
  
  # Get the function name
  fun_name <- as.character(fun_call[[1]])
  vebcat(fun_name, color = "debugLight")
  
  # Get the arguments
  args_list <- as.list(fun_call)[-1]
  
  # Get the environment of the function that called this function
  fun_env <- parent.frame(2)
  
  for (arg in names(args_list)) {
    catn(arg)
    print(class(arg))
    print(head(eval(args_list[[arg]], envir = fun_env), 3))
  }
  
  catn()
}

vebprint <- function(x, veb = TRUE, text = NULL) {
  if (veb) {
    if (!is.null(text)) {
      cat(text, "\n")
    }
    print(x)
  }
}

catn <- function(...) {
  args <- list(...)
  text <- do.call(paste, c(list(sep = " "), as.character(args)))
  
  cat(text, "\n")
}

colcat <- function(..., color = NULL, veb = TRUE) {
  if (veb) {
    args <- list(...)
    text <- do.call(paste, c(list(sep = " "), as.character(args)))
    
    if (!is.null(color)) {
      return(paste(cc[[color]](text)))
    } else {
      return(text)
    }
  }
}

vebcat <- function(..., color = NULL, veb = TRUE) {
  if (veb) {
    args <- list(...)
    text <- do.call(paste, c(list(sep = " "), as.character(args)))
    
    if (!is.null(color)) {
      cat(cc[[color]](text), "\n")
    } else {
      cat(text, "\n")
    }
  }
}

highcat <- function(...) {
  args <- list(...)
  text <- do.call(paste, c(list(sep = " "), as.character(args)))
  
  return(cc$highlight(text))
}
