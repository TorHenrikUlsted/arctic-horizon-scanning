message("Initiating AMBIO wrangling")
# Here we wrangle with the AMBIO data to make it more R friendly
## Read file
ambio_preformatted = read.csv("./resources/ambio.csv", header = F)

# conduct the wrangling

## remove the not important rows and columns
ambio_formatted = ambio_preformatted[-c(1, 3:9), ]
ambio_formatted = ambio_formatted[,-25]

## set the first row to be the header and remove the row from the dataset
colnames(ambio_formatted) = as.character(unlist(ambio_formatted[1, ]))
ambio_formatted = ambio_formatted[-1, ]
## give the first column a name
colnames(ambio_formatted)[1] = "Species_SubSpecies"


# Create present and absent lists
## create conditions
### symbols meaning: Present: ●, IR, IT and Absent: ○, ?, × 
condition1 = apply(ambio_formatted[, 2:24], 1, function(x) any(x %in% c("●", "IR", "IT")))
## Use the condition to create present and absent lists
ambio_present = subset(ambio_formatted, subset = condition1)
ambio_present = select(ambio_present, Species_SubSpecies)
ambio_present$Species_SubSpecies = trimws(ambio_present$Species_SubSpecies)

## Only outputs unquie species names
ambio_absent = subset(ambio_formatted, subset = !condition1)
ambio_absent = select(ambio_absent, Species_SubSpecies)
ambio_absent$Species_SubSpecies = trimws(ambio_absent$Species_SubSpecies)

write.csv(ambio_preformatted, "./outputs/filtering/origin_lists/ambio_preformatted.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(ambio_formatted, "./outputs/filtering/wrangling_outputs/ambio/ambio_formatted.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(ambio_present, "./outputs/filtering/wrangling_outputs/ambio/ambio_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(ambio_absent, "./outputs/filtering/wrangling_outputs/ambio/ambio_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")

cat("AMBIO data wrangling completed \n")