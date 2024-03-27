wrangle_aba <- function(column, verbose = F) {
  vebcat("Initiating ABA wrangling protocol", color = "funInit")
  # read CSV file
  aba_preformat <- fread("./resources/data-raw/aba.csv", header = F)

  ## format the ABA CSV file
  # Remove empty columns
  vebcat("selecting columns", veb = verbose)
  aba_selected <- select(aba_preformat, -tail(seq_along(aba_preformat), 3))

  # Assign new column names using two rows
  colnames(aba_selected) <- paste(aba_selected[3, ], aba_selected[5, ])
  # Trim spaces in the headings
  colnames(aba_selected) <- trimws(colnames(aba_selected))
  # Add new names to columns 29 to 42
  colnames(aba_selected)[29:39] <- c("arcticOccurence", "arcticEndemicSpeciesAE", "borderline", "introduced", "naturalized", "nonNativeStableCasual", "stableCasual", "nativeCasual", "paf", "genusCount", "familyCount")
  # Remove row 1:6
  aba_selected <- aba_selected[-c(1:6), ]

  # Assume that aba_selected is your existing data frame
  # and that the first column contains the mixed Family, Genus, and Species text

  # Create new columns for Class, Family, Genus, and Species
  aba_selected$class <- ""
  aba_selected$family <- ""
  aba_selected$genus <- ""
  aba_selected$species <- ""
  aba_selected$subSpecies <- ""

  # Initialize variables to keep track of the current classification levels
  current_class <- ""
  current_family <- ""
  current_genus <- ""
  
  vebcat("Sorting names.", veb = verbose)
    
  for (i in seq_len(nrow(aba_selected))) {
    cat("\rSequencing row:", cc$lightSteelBlue(i), " / ", highcat(nrow(aba_selected)), fill = F)
    flush.console()
    # Get the text from the first column
    line <- as.character(aba_selected[i, 1])

    # Split the line into its components
    components <- strsplit(line, " ")[[1]]
    # Check the number of components to determine the classification level
    if (length(components) == 1) {
      # Class level
      current_class <- components[1]
      current_family <- ""
      current_genus <- ""

      # skip the first row of every class as it has no info other than the class name
      if (aba_selected[i, "class"] != "") {
        aba_selected[i, "class"] <- current_class
      }
    } else if (length(components) >= 2) {
      # Family, Genus, or Species level

      # Check the number of digits in the first component
      num_digits <- nchar(gsub("[^0-9]", "", components[1]))

      if (num_digits == 2) {
        # Family level
        current_family <- components[2]
        current_genus <- ""

        # Update the values in the data frame
        aba_selected[i, "class"] <- current_class
        aba_selected[i, "family"] <- current_family
      } else if (num_digits == 4) {
        # Genus level
        current_genus <- components[2]

        # Update the values in the data frame
        aba_selected[i, "class"] <- current_class
        aba_selected[i, "family"] <- current_family
        aba_selected[i, "genus"] <- current_genus
      } else if (num_digits >= 6) {
        # Species level

        # Remove "ssp." prefix and any first capital letter followed by a dot from species name
        species_name <- gsub("^(ssp\\.\\s)?([A-Z]\\.\\s)?", "", paste(components[2:length(components)], collapse = " "))

        # Check if the species name contains a sub-species designation with or without parenthesis
        if (grepl("\\(ssp\\.\\s.*\\)", species_name)) {
          # If so, assign the sub-species to the subSpecies column
          sub_species <- gsub("^.*\\(ssp\\.\\s(.*)\\)$", "\\1", species_name)
          aba_selected[i, "subSpecies"] <- sub_species

          # Remove the sub-species designation from the species name
          species_name <- gsub("\\(ssp\\.\\s.*\\)", "", species_name)
        } else if (grepl("ssp\\.\\s", species_name)) {
          # If so, assign the sub-species to the subSpecies column
          sub_species <- gsub("^(.*)(ssp\\.\\s)(.*)$", "\\3", species_name)
          aba_selected[i, "subSpecies"] <- sub_species

          # Remove the sub-species designation from the species name
          species_name <- gsub("^(.*)(\\s+ssp\\.\\s+.*)$", "\\1", species_name)
        }

        # Update the values in the data frame
        aba_selected[i, "class"] <- current_class
        aba_selected[i, "family"] <- current_family
        aba_selected[i, "genus"] <- current_genus
        aba_selected[i, "species"] <- species_name
      }
    }
  };catn()
  
  vebcat("Formatting df.", veb = verbose)
  ## Remove column 1 which is now the old information of class, family, genus, species and subSpecies
  aba_formatted <- aba_selected[, -1]
  ## Remove empty rows and rows with all columns NA
  aba_formatted <- aba_formatted[!apply(aba_formatted, 1, function(x) all(is.na(x) | x == "")), ]
  ## Create new column with Genus and species Combined
  aba_formatted$scientificName <- with(aba_formatted, paste(genus, species, subSpecies, sep = " "))
  names(aba_formatted)
  ## Move the new columns to be the first four columns of the data frame
  col_names <- colnames(aba_formatted)
  new_order <- c(col_names[(length(col_names) - 5):length(col_names)], col_names[1:(length(col_names) - 6)])
  aba_formatted <- aba_formatted[, ..new_order]
  ## Change to lowercase in class and family for cleaner look
  aba_formatted$class <- toupper(substr(aba_formatted$class, 1, 1)) %>% paste0(tolower(substr(aba_formatted$class, 2, nchar(aba_formatted$class))))
  aba_formatted$family <- toupper(substr(aba_formatted$family, 1, 1)) %>% paste0(tolower(substr(aba_formatted$family, 2, nchar(aba_formatted$family))))
  ## Replace the square symbol with dash in all columns
  aba_formatted[] <- lapply(aba_formatted, function(x) gsub("", "-", x))

  # Distinguish between present and absent data
  ## create a new dataset with species absent from the Arctic
  
  vebcat("Creating conditions.", veb = verbose)
  ## Create conditions
  c1 <- ifelse(aba_formatted$borderline == 1, TRUE, FALSE)
  c2 <- apply(aba_formatted[, 7:32], 1, function(x) all(x %in% c("-", "?", "**")))
  
  vebcat("Making ABA present.", veb = verbose)
  ## use the !conditions to get species present in the Arctic
  aba_present <- aba_formatted[!(c1 | c2), ]
  ## remove the rows with empty cell from column 7 to 33.
  aba_present <- aba_present[!apply(aba_present[, 7:33], 1, function(x) any(x == "")), ]
  ## Select only the species_SubSpecies column
  aba_present <- select(aba_present, scientificName)
  aba_present$scientificName <- trimws(aba_present$scientificName)
  
  setnames(aba_present, old = "scientificName", new = column)

  vebcat("Making ABA absent.", veb = verbose)
  ## use conditions to create a new dataset by choosing only those satisfying the conditions from the aba_formatted dataset
  aba_absent <- aba_formatted[c1 | c2, ]
  ## remove the rows with empty cell from column 7 to 33.
  aba_absent <- aba_absent[!apply(aba_absent[, 7:33], 1, function(x) any(x == "")), ]
  ## select only the species with subspecies name
  aba_absent <- select(aba_absent, scientificName)
  aba_absent$scientificName <- trimws(aba_absent$scientificName)
  
  setnames(aba_absent, old = "scientificName", new = column)

  vebcat("Writing out files.", veb = verbose)
  
  create_dir_if("./outputs/setup/wrangle/aba")
  
  set_df_utf8(aba_formatted)
  fwrite(aba_formatted, "./outputs/setup/wrangle/aba/aba-formatted.csv", row.names = F, bom = T)

  set_df_utf8(aba_present)
  fwrite(aba_present, "./outputs/setup/wrangle/aba/aba-present.csv", row.names = F, bom = T)
  
  set_df_utf8(aba_absent)
  fwrite(aba_absent, "./outputs/setup/wrangle/aba/aba-absent.csv", row.names = F, bom = T)

  vebcat("ABA wrangling protocol completed successfully.", color = "funSuccess")
  
  return(list(
    aba_present = aba_present,
    aba_absent = aba_absent
  ))
}
