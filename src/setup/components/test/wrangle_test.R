wrangle_test_small <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/test")
  create_dir_if(dir)
  
  absent_out <- paste0(dir, "/", name, "/", name, "-absent.csv")
  
  formatted <- fread(paste0("./resources/data-raw/test/", name, "/", name, ".csv"), sep = "\t")
  
  setnames(formatted, old = "scientificName", new = "verbatimName")
  
  formatted[, (column) := verbatimName]
  
  absent <- formatted
  
  fwrite(absent, absent_out, bom = T)
  
  write_wrangled_md(
    dt.list = list(formatted = formatted, absent = absent),
    name = "Test Small",
    column = column
  )
  
  return(list(
    absent = absent
  ))
}

wrangle_test_big <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/test")
  create_dir_if(dir)
  
  absent_out <- paste0(dir, "/", name, "/", name, "-absent.csv")
  
  formatted <- fread(paste0("./resources/data-raw/test/", name, "/", name, ".csv"))
  
  setnames(formatted, old = "scientificName", new = "verbatimName")
  setnames(formatted, old = "author", new = "verbatimNameAuthorship")
  
  formatted[, (column) := verbatimName]
  formatted[, (paste0(column, "Authorship")) := verbatimNameAuthorship]
  absent <- formatted
  
  fwrite(absent, absent_out, bom = T)
  
  write_wrangled_md(
    dt.list = list(formatted = formatted, absent = absent),
    name = "Test Big",
    column = column
  )
  
  return(list(
    absent = absent
  ))
}

wrangle_test_known <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/test")
  create_dir_if(dir)
  
  present_out <- paste0(dir, "/", name, "/", name, "-present.csv")
  
  formatted <- fread(paste0("./resources/data-raw/test/", name, "/", name, ".csv"), sep = "\t")
  
  setnames(formatted, old = "scientificName", new = "verbatimName")
  
  formatted[, (column) := verbatimName]
  
  present <- formatted
  
  fwrite(present, present_out, bom = T)
  
  write_wrangled_md(
    dt.list = list(formatted = formatted, present = present),
    name = "Test known",
    column = column
  )
  
  return(list(
    present = present
  ))
}
