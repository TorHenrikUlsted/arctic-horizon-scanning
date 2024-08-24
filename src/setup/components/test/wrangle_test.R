wrangle_test_small <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/test")
  create_dir_if(dir)
  
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
  formatted <- fread(paste0("./resources/data-raw/test/", name, ".csv"), sep = "\t")
  
  setnames(formatted, old = "scientificName", new = column)
  
  absent <- set_df_utf8(formatted)
  fwrite(absent, absent_out, row.names = F, bom = T)
  
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
  
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
  formatted <- fread(paste0("./resources/data-raw/test/", name, ".csv"), sep = "\t")
  
  setnames(formatted, old = "scientificName", new = column)
  
  absent <- set_df_utf8(formatted)
  fwrite(absent, absent_out, row.names = F, bom = T)
  
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
  
  present_out <- paste0(dir, "/", name, "-present.csv")
  
  formatted <- fread(paste0("./resources/data-raw/test/", name, ".csv"), sep = "\t")
  
  setnames(formatted, old = "scientificName", new = column)
  
  present <- set_df_utf8(formatted)
  fwrite(present, present_out, row.names = F, bom = T)
  
  write_wrangled_md(
    dt.list = list(formatted = formatted, present = present),
    name = "Test known",
    column = column
  )
  
  return(list(
    present = present
  ))
}
