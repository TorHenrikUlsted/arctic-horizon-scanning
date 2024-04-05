# wrangle_template is ignored, change the name
wrangle_template <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)
  
  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  present_out <- paste0(dir, "/", name, "-present.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
  dt <- fread(paste0("./resources/data-raw/", name, ".csv"))
  
  ## ADD wrangling in here
  
  return(list(
    present = present, # Remove if you only have one return
    absnt = absent
  ))
}