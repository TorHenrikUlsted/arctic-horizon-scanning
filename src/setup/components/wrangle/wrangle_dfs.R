wrangle_dfs <- function(src, column) {
  files <- list.files(path = src, pattern = "\\.R$")

  # Source all files
  lapply(paste0(src, files), source)

  results <- list()

  # Call their functions
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
        cat(green(toString(class(get(name))), "\n"))
      }
    }

    results[[name]] <- get(name)
  }

  return(results)
}
