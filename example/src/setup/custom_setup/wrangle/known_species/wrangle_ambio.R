wrangle_ambio <- function(name, column, verbose = FALSE) {
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)

  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  present_out <- paste0(dir, "/", name, "-present.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")

  preformatted <- fread(paste0("./resources/data-raw/", name, ".csv"), header = F)

  vebcat("filtering rows.", veb = verbose)
  ## remove the not important rows and columns
  formatted <- preformatted[-c(1, 3:9), ]
  formatted <- formatted[, -25]

  ## set the first row to be the header and remove the row from the dataset
  colnames(formatted) <- as.character(unlist(formatted[1, ]))
  formatted <- formatted[-1, ]
  ## give the first column a name
  colnames(formatted)[1] <- "scientificName"
  formatted <- unique(formatted, by = "scientificName")

  vebcat("Creating present df.", veb = verbose)
  # Create present and absent lists
  ## create conditions
  ### symbols meaning: Present: ●, IR, IT and Absent: ○, ?, ×
  condition1 <- apply(formatted[, 2:24], 1, function(x) any(x %in% c("●", "IR", "IT")))
  ## Use the condition to create present and absent lists
  present <- subset(formatted, subset = condition1)
  present <- select(present, scientificName)
  present$scientificName <- trimws(present$scientificName)

  setnames(present, old = "scientificName", new = column)

  vebcat("Creating absent df.", veb = verbose)
  ## Only outputs unquie species names
  absent <- subset(formatted, subset = !condition1)
  absent <- select(absent, scientificName)
  absent$scientificName <- trimws(absent$scientificName)

  setnames(absent, old = "scientificName", new = column)

  catn("Writing out files to:", colcat(dir, color = "output"))

  formatted <- set_df_utf8(formatted)
  fwrite(formatted, formatted_out, row.names = F, bom = T)

  present <- set_df_utf8(present)
  fwrite(present, present_out, row.names = F, bom = T)

  absent <- set_df_utf8(absent)
  fwrite(absent, absent_out, row.names = F, bom = T)
  
  fl <- nrow(formatted)
  pl <- nrow(present)
  al <- nrow(absent)
  
  md_dt <- data.table(
    formatted = fl,
    present = pl,
    absent = al,
    lost = fl - sum(pl, al)
  )
  
  mdwrite(
    config$files$post_seq_md,
    text = "3;Wrangled AMBIO",
    data = md_dt
  )
  
  return(list(
    present = present,
    absent = absent
  ))
}
