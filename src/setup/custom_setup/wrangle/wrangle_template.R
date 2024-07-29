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
  
  
  ## ADD MD
  mdwrite(
    post_seq_nums,
    text = paste0(
      "2;Template\n\n",
      "Number of species in Template formatted:", nrow(formatted), "**  ",
      "Number of species in Template present:", nrow(template_present), "**   ",
      "Number of species in Template absent:", nrow(template_absent), "**   ",
    )
  )
  
  return(list(
    present = present, # Remove if you only have one return
    absnt = absent
  ))
}