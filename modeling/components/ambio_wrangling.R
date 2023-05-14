# Here we wrangle with the AMBIO data to make it more R friendly
## Read file
ambio_preformatted = read.csv("resources/AMBIO.csv", header = F)

# conduct the wrangling

## remove the not important rows and columns
ambio_formatted = ambio_preformatted[-c(1, 3:10), ]
ambio_formatted = ambio_formatted[,-25]

## set the first row to be the header and remove the row from the dataset
colnames(ambio_formatted) = as.character(unlist(ambio_formatted[1, ]))
ambio_formatted = ambio_formatted[-1, ]
## give the first column a name
colnames(ambio_formatted)[1] = "Species_SubSpecies"


# Create present and absent lists
## create conditions
### symbols meaning: Present: ●, I, IT and Absent: ○, ?, × 
condition1 <- apply(ambio_formatted[, 2:24], 1, function(x) any(x %in% c("●", "I", "IT")))
## Use the condition to create present and absent lists
ambio_arctic_present <- subset(ambio_formatted, subset = condition1)
ambio_arctic_absent <- subset(ambio_formatted, subset = !condition1)

## Only outputs unquie species names
ambio_arctic_present = select(ambio_arctic_present, Species_SubSpecies)
ambio_arctic_absent = select(ambio_arctic_absent, Species_SubSpecies)

# Conduct WFO synonym matching

## Check if response is "y"
if (abaAmbioList_response == "y") {
  ### Do not run this part of the script
  cat("Skipping the WFO synonym check for the AMBIO list. \n")
} else {
  ### Run this part of the script
  cat("Running the WFO synonym check on the AMBIO list.\n")
  ## Synonym check for the Arctic present species
  ### Take the time
  start_time = Sys.time()
  wfo_ambio_arctic_present = WFO.match(spec.data = ambio_arctic_present, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 500)
  end_time = Sys.time()
  ### Calculate the time
  formatted_elapsed_time = format_elapsed_time(start_time, end_time)
  ## Print message with elapsed time
  cat("WFO completed the match for ambio_present species in ", formatted_elapsed_time, "\n")
  ## write it into a CSV file
  cat("Creating CSV file of the output \n")
  write.csv(wfo_ambio_arctic_present, "outputs/wfo_ambio_arctic_present.csv")
  ## Choose the scientifically accepted name 
  ambio_arctic_present = select(wfo_ambio_arctic_present, scientificName)
  ## Remove duplicates
  ambio_arctic_present = distinct(ambio_arctic_present)
  
  ## Synonym check for Arctic absent species
  ### Take the time
  start_time = Sys.time()
  wfo_ambio_arctic_absent = WFO.match(spec.data = ambio_arctic_absent, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 500)
  end_time = Sys.time()
  ### Calculate the time
  formatted_elapsed_time = format_elapsed_time(start_time, end_time)
  ## Print message with elapsed time
  cat("WFO completed the match for ambio_absent species in ", formatted_elapsed_time, "\n")
  ## write it into a CSV file
  cat("Creating CSV file of the output \n")
  write.csv(wfo_ambio_arctic_absent, "outputs/wfo_ambio_arctic_absent.csv")
  ## Choose the scientifically accepted name 
  ambio_arctic_absent = select(wfo_ambio_arctic_absent, scientificName)
  ## Remove duplicates
  ambio_arctic_absent = distinct(ambio_arctic_absent)
}

cat("AMBIO data wrangling complete \n")