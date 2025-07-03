# 'wrangle_template' is ignored, change the name to 'wrangle_newname'
wrangle_template <- function(name, column, verbose = FALSE) {
  # The name is taken from the function name
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)

  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  present_out <- paste0(dir, "/", name, "-present.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")

  if (file.exists(formatted_out) && file.exists(present_out) && file.exists(absent_out)) {
    present <- fread(present_out)
    absent <- fread(absent_out)
  } else {
    dt <- fread(paste0("./resources/data-raw/", name, ".csv"))
    ## ADD wrangling in here
    
    ## name the final result objects present & absent
    
    
    
    ## ADD MD
    mdwrite(
      config$files$post_seq_md,
      text = paste0(
        "2;Template\n\n",
        "Number of species in Template formatted:", nrow(formatted), "**  ",
        "Number of species in Template present:", nrow(present), "**   ",
        "Number of species in Template absent:", nrow(absent), "**   ",
      )
    )
  } # End if files exist else statement
  
  present[, sourceDataset := name]
  absent[, sourceDataset := name]

  return(list(
    present = present,
    absent = absent
  ))
}
