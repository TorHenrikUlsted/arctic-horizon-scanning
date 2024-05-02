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

warn <- function(w, warn.file, warn.txt, iteration = NULL) {
  warn_msg <- conditionMessage(w)
  create_file_if(warn.file, keep = TRUE)
  warn_con <- file(warn.file, open = "a")
  
  if (!is.null(iteration)) {
    writeLines(paste(warn.txt, "in iteration", iteration, ":", warn_msg), warn_con)
  } else {
    writeLines(paste(warn.txt, ":", warn_msg), warn_con)
  }
  
  close(warn_con)
  invokeRestart(findRestart("muffleWarning"))
}

err <- function(e, err.file, err.txt, iteration = NULL) {
  err_msg <- conditionMessage(e)
  create_file_if(err.file, keep = TRUE)
  err_con <- file(err.file, open = "a")
  
  if (!is.null(iteration)) {
    writeLines(paste(err.txt, "in iteration", iteration, ":", err_msg), err_con)
  } else {
    writeLines(paste(err.txt, ":", err_msg), err_con)
  }
  close(err_con)
  stop(e)
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
      if (color == "funInit" || color == "seqInit" || color == "proInit") {
        cat("\n", cc[[color]](text), "\n")
      } else if (color == "funSuccess" || color == "seqSuccess" || color == "proSuccess") {
        cat("", cc[[color]](text), "\n\n")
      } else {
        cat(cc[[color]](text), "\n")  
      }
      
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
