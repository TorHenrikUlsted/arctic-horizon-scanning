wrangle_ambio <- function(column, verbose = F) {
  vebcat("Initiating AMBIO wrangling protocol.", color = "funInit")

  ambio_preformatted <- read.csv("./resources/data-raw/ambio.csv", header = F)
  
  vebcat("filtering rows.", veb = verbose)
  ## remove the not important rows and columns
  ambio_formatted <- ambio_preformatted[-c(1, 3:9), ]
  ambio_formatted <- ambio_formatted[, -25]

  ## set the first row to be the header and remove the row from the dataset
  colnames(ambio_formatted) <- as.character(unlist(ambio_formatted[1, ]))
  ambio_formatted <- ambio_formatted[-1, ]
  ## give the first column a name
  colnames(ambio_formatted)[1] <- "scientificName"

  vebcat("Creating present df.", veb = verbose)
  # Create present and absent lists
  ## create conditions
  ### symbols meaning: Present: ●, IR, IT and Absent: ○, ?, ×
  condition1 <- apply(ambio_formatted[, 2:24], 1, function(x) any(x %in% c("●", "IR", "IT")))
  ## Use the condition to create present and absent lists
  ambio_present <- subset(ambio_formatted, subset = condition1)
  ambio_present <- select(ambio_present, scientificName)
  ambio_present$scientificName <- trimws(ambio_present$scientificName)
  
  setnames(ambio_present, old = "scientificName", new = column)

  vebcat("Creating absent df.", veb = verbose)
  ## Only outputs unquie species names
  ambio_absent <- subset(ambio_formatted, subset = !condition1)
  ambio_absent <- select(ambio_absent, scientificName)
  ambio_absent$scientificName <- trimws(ambio_absent$scientificName)
  
  setnames(ambio_absent, old = "scientificName", new = column)
  
  vebcat("Writing out files.", veb = verbose)
  create_dir_if("./outputs/setup/wrangle/ambio")
  
  ambio_formatted <- set_df_utf8(ambio_formatted)
  fwrite(ambio_formatted, "./outputs/setup/wrangle/ambio/ambio-formatted.csv", row.names = F, bom = T)
  
  ambio_present <- set_df_utf8(ambio_present)
  fwrite(ambio_present, "./outputs/setup/wrangle/ambio/ambio-present.csv", row.names = F, bom = T)
  
  ambio_absent <- set_df_utf8(ambio_absent)
  fwrite(ambio_absent, "./outputs/setup/wrangle/ambio/ambio-absent.csv", row.names = F, bom = T)

  vebcat("AMBIO wrangling protocol successfully completed.", color = "funSuccess")

  return(list(
    ambio_present = ambio_present,
    ambio_absent = ambio_absent
  ))
}
