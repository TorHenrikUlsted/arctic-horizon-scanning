wrangle_test_small <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)
  
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
  dt <- fread(paste0("./resources/data-raw/test/", name, ".csv"), sep = "\t")
  
  setnames(dt, old = "scientificName", new = column)
  
  dt <- set_df_utf8(dt)
  fwrite(dt, absent_out, row.names = F, bom = T)
  
  return(list(
    absent = dt
  ))
}

wrangle_test_big <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)
  
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
  dt <- fread(paste0("./resources/data-raw/test/", name, ".csv"), sep = "\t")
  
  setnames(dt, old = "scientificName", new = column)
  
  dt <- set_df_utf8(dt)
  fwrite(dt, absent_out, row.names = F, bom = T)
  
  return(list(
    absent = dt
  ))
}

wrangle_test <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)
  
  present_out <- paste0(dir, "/", name, "-present.csv")
  
  dt <- fread(paste0("./resources/data-raw/test/", name, ".csv"), sep = "\t")
  
  setnames(dt, old = "scientificName", new = column)
  
  dt <- set_df_utf8(dt)
  fwrite(dt, present_out, row.names = F, bom = T)
  
  return(list(
    present = dt
  ))
}
