wrangle_dfs <- function(src, column, dynamic.name.source) {
  files <- list.files(path = src, pattern = "\\.R$")

  # Source all files
  lapply(paste0(src, files), source)

  results <- list()
  
  if (dynamic.name.source == "file") {
    for (file in files) {
      name <- sub("\\.R$", "", file)
      
      # create a name for each
      assign(name, NULL)
      
      func_name <- paste0("wrangle_", name)
      
      if (exists(func_name)) {
        assign(name, do.call(func_name, list(column = column)))
        
        # Check if the result is a data.table or data.frame
        if (!"data.table" %in% class(get(name)) && !"data.frame" %in% class(get(name))) {
          stop(paste("The result of", func_name, "is not a 'data.table' or 'data.frame'."))
        } else {
          catn(toString(class(get(name))))
        }
      }
      
      results[[name]] <- get(name)
    }
    
  } else if (dynamic.name.source == "object") {
    for (name in names(results)) {
      # Check if the object exists
      if (exists(name)) {
        # Get the object
        obj <- get(name)
        
        # Check if the object is a data.table or data.frame
        if (!"data.table" %in% class(obj) && !"data.frame" %in% class(obj)) {
          stop(paste("The object", name, "is not a 'data.table' or 'data.frame'."))
        } else {
          catn(paste("The class of the object", name, "is", class(obj)))
        }
        
        # Assign the object to the results list
        results[[name]] <- obj
      }
    }
  }

  return(results)
}
