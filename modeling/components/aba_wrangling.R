# read CSV file
ABA_preformat = read.csv("resources/ABA_2013.csv", header = F)
ncol(ABA_preformat)

## format the ABA CSV file
# Remove empty columns
ABA_formatted = select(ABA_preformat, -tail(seq_along(ABA_preformat), 3))

# Assign new column names using two rows
colnames(ABA_formatted) = paste(ABA_formatted[3, ], ABA_formatted[5, ])
# Add new names to columns 29 to 42
colnames(ABA_formatted)[29:39] = c("ArcticOccurence", "ArcticEndemicSpeciesAE", "Borderline", "Introduced", "Naturalized", "nonNativeStableCasual", "StableCasual", "NativeCasual", "PAF", "GenusCount", "FamilyCount")
ABA_formatted = ABA_formatted[-c(3,5), ]
# Remove row 1:6
ABA_formatted = ABA_formatted[-c(1:4), ]

# Assume that ABA_formatted is your existing data frame
# and that the first column contains the mixed Family, Genus, and Species text

# Create new columns for Class, Family, Genus, and Species
ABA_formatted$Class = ""
ABA_formatted$Family = ""
ABA_formatted$Genus = ""
ABA_formatted$Species = ""
ABA_formatted$Subspecies = ""

# Initialize variables to keep track of the current classification levels
current_class = ""
current_family = ""
current_genus = ""

for (i in seq_len(nrow(ABA_formatted))) {
  
  # Get the text from the first column
  line = ABA_formatted[i,1]
  
  # Split the line into its components
  components = strsplit(line, " ")[[1]]
  
  # Check the number of components to determine the classification level
  if (length(components) == 1) {
    # Class level
    current_class = components[1]
    current_family = ""
    current_genus = ""
    
    #skip the first row of every class as it has no info other than the class name
    if (ABA_formatted[i,"Class"] != "") {
      ABA_formatted[i,"Class"] = current_class
    }
    
  } else if (length(components) >= 2) {
    # Family, Genus, or Species level
    
    # Check the number of digits in the first component
    num_digits = nchar(gsub("[^0-9]", "", components[1]))
    
    if (num_digits == 2) {
      # Family level
      current_family = components[2]
      current_genus = ""
      
      # Update the values in the data frame
      ABA_formatted[i,"Class"] = current_class
      ABA_formatted[i,"Family"] = current_family
      
    } else if (num_digits == 4) {
      # Genus level
      current_genus = components[2]
      
      # Update the values in the data frame
      ABA_formatted[i,"Class"] = current_class
      ABA_formatted[i,"Family"] = current_family
      ABA_formatted[i,"Genus"] = current_genus
    } else if (num_digits >= 6) {
      # Species level
      
      # Remove "ssp." prefix and any first capital letter followed by a dot from species name
      species_name = gsub("^(ssp\\.\\s)?([A-Z]\\.\\s)?", "", paste(components[2:length(components)], collapse = " "))
      
      # Check if the species name contains a sub-species designation with or without parenthesis
      if (grepl("\\(ssp\\.\\s.*\\)", species_name)) {
        # If so, assign the sub-species to the subSpecies column
        sub_species = gsub("^.*\\(ssp\\.\\s(.*)\\)$", "\\1", species_name)
        ABA_formatted[i,"Subspecies"] = sub_species
        
        # Remove the sub-species designation from the species name
        species_name = gsub("\\(ssp\\.\\s.*\\)", "", species_name)
      } else if (grepl("ssp\\.\\s", species_name)) {
        # If so, assign the sub-species to the subSpecies column
        sub_species = gsub("^(.*)(ssp\\.\\s)(.*)$", "\\3", species_name)
        ABA_formatted[i,"Subspecies"] = sub_species
        
        # Remove the sub-species designation from the species name
        species_name = gsub("^(.*)(\\s+ssp\\.\\s+.*)$", "\\1", species_name)
      }
      
      # Update the values in the data frame
      ABA_formatted[i,"Class"] = current_class
      ABA_formatted[i,"Family"] = current_family
      ABA_formatted[i,"Genus"] = current_genus
      ABA_formatted[i,"Species"] = species_name
    }
  }
}

## Remove column 1 which is now the old information of class, family, genus, species and subSpecies
ABA_formatted = ABA_formatted[, -1]
## Remove empty rows and rows with all columns NA
ABA_formatted = ABA_formatted[!apply(ABA_formatted, 1, function(x) all(is.na(x) | x == "")),]
## Create new column with Genus and species Combined
ABA_formatted$Species_SubSpecies = with(ABA_formatted, paste (Genus, Species, Subspecies, sep = " "))
## Move the new columns to be the first four columns of the data frame
ABA_formatted = ABA_formatted[,c("Class", "Family", "Genus", "Species", "Subspecies", "Species_SubSpecies", setdiff(colnames(ABA_formatted), c("Class", "Family", "Genus", "Species", "Subspecies", "Species_SubSpecies")))]
## Change to lowercase in class and family for cleaner look
ABA_formatted$Class = toupper(substr(ABA_formatted$Class, 1, 1)) %>% paste0(tolower(substr(ABA_formatted$Class, 2, nchar(ABA_formatted$Class))))
ABA_formatted$Family = toupper(substr(ABA_formatted$Family, 1, 1)) %>% paste0(tolower(substr(ABA_formatted$Family, 2, nchar(ABA_formatted$Family))))
## Replace the square symbol with dash in all columns
ABA_formatted[] = lapply(ABA_formatted, function(x) gsub("", "-", x))

# Distinguish between present and absent data
## create a new dataset with species absent from the Arctic

## Create conditions
c1 = ifelse(ABA_formatted$Borderline == 1, TRUE, FALSE)
c2 = apply(ABA_formatted[, 7:32], 1, function(x) all(x %in% c("-", "?", "**")))

## use the !conditions to get species present in the Arctic
aba_arctic_present = ABA_formatted[!(c1 | c2), ]
## remove the rows with empty cell from column 7 to 33.
aba_arctic_present = aba_arctic_present[!apply(aba_arctic_present[, 7:33], 1, function(x) any(x == "")), ]
## Select only the species_SubSpecies column
aba_arctic_present = select(aba_arctic_present, Species_SubSpecies)

## use conditions to create a new dataset by choosing only those satisfying the conditions from the ABA_formatted dataset
aba_arctic_absent = ABA_formatted[c1 | c2, ]
## remove the rows with empty cell from column 7 to 33.
aba_arctic_absent = aba_arctic_absent[!apply(aba_arctic_absent[, 7:33], 1, function(x) any(x == "")), ]
## select only the species with subspecies name
aba_arctic_absent = select(aba_arctic_absent, Species_SubSpecies)

# Check if response is "y"
if (abaAmbioList_response == "y") {
  # Do not run this part of the script
  cat("Skipping the WFO synonym check for the ABA list. \n")
} else {
  # Run this part of the script
  cat("Running the WFO synonym check for the ABA list.\n")
  
  # Conduct WFO synonym matching
  ## Synonym check for the Arctic present species
  ### take the time
  start_time = Sys.time()
  wfo_aba_arctic_present = WFO.match(spec.data = aba_arctic_present, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 500)
  end_time = Sys.time()
  ### Calculate the time
  formatted_elapsed_time = format_elapsed_time(start_time, end_time)
  ## Print message with elapsed time
  cat("WFO completed the match for aba_present species in ", formatted_elapsed_time, "\n")
  ## write it into a CSV file
  cat("Creating CSV file of the output \n")
  write.csv(wfo_aba_arctic_present, "outputs/wfo_aba_arctic_present.csv")
  ## Choose the scientifically accepted name 
  aba_arctic_present = select(wfo_aba_arctic_present, scientificName)
  ## Remove duplicates
  aba_arctic_present = distinct(aba_arctic_present)
  
  ## Synonym check for Arctic absent species
  ### Take time
  start_time = Sys.time()
  wfo_aba_arctic_absent = WFO.match(spec.data = aba_arctic_absent, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 500)
  end_time = Sys.time()
  ### Calculate the time
  formatted_elapsed_time = format_elapsed_time(start_time, end_time)
  ## Print message with elapsed time
  cat("WFO completed the match for aba_absent species in ", formatted_elapsed_time, "\n")
  ## write it into a CSV file
  cat("Creating CSV file of the output \n")
  write.csv(wfo_aba_arctic_absent, "outputs/wfo_aba_arctic_absent.csv")
  ## Choose the scientifically accepted name 
  aba_arctic_absent = select(wfo_aba_arctic_absent, scientificName)
  ## Remove duplicates
  aba_arctic_absent = distinct(aba_arctic_absent)
}

cat("ABA data wrangling complete \n")