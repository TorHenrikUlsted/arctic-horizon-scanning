wrangle_aba <- function(name, column, verbose = F) {
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)

  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  present_out <- paste0(dir, "/", name, "-present.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")

  # read CSV file
  
  preformat <- fread(paste0("./resources/data-raw/", name, ".csv"), header = F)
  
  ## format the ABA CSV file
  vebcat("selecting columns", veb = verbose)
  aba_selected <- preformat
  # Assign new column names using two rows
  data.table::setnames(aba_selected, paste(aba_selected[3], aba_selected[5]))
  # Trim spaces in the headings
  setnames(aba_selected, trimws(names(aba_selected)))
  
  # Add new names to columns 29 to 39
  new_names <- c("arcticOccurence", "arcticEndemicSpeciesAE", "borderline", "introduced", "naturalized", "nonNativeStableCasual", "stableCasual", "nativeCasual", "paf", "genusCount", "familyCount")
  setnames(aba_selected, 29:39, new_names)
  
  # Remove the last three NA columns
  aba_selected <- aba_selected[, -((ncol(aba_selected) - 2):ncol(aba_selected)), with = FALSE]

  # Remove row 1:6
  aba_selected <- aba_selected[7:.N]
  setnames(aba_selected, 1, "verbatimName")
  aba_selected[, speciesCode := verbatimName]
  # Remove text in the end and blank rows
  aba_selected <- aba_selected[!grepl("^Total|^Number|^Mean", speciesCode) & speciesCode != ""]
  # Remove parenthesis in the
  aba_selected <- aba_selected[, speciesCode := gsub("\\(|\\)", "", speciesCode)]
  
  aba_selected[, speciesCode := {
    tmp <- speciesCode
    tmp <- clean_string(tmp, verbose)
    tmp <- clean_designations(tmp, config$species$standard_infraEpithets, verbose)
    clean_symbols(tmp, config$species$standard_symbols, verbose)
  }]
  
  vebcat("Sorting names.", veb = verbose)
  
  # Create new columns for Class, Family, Genus, and Species
  preformat[, `:=`(
    class = "",
    family = "",
    genus = "",
    species = "",
    infraEpithet = "",
    scientificName = ""
  )]
  
  # Initialize variables to keep track of current taxonomic levels
  current_class <- ""
  current_family <- ""
  current_genus <- ""
  
  process_row <- function(col, current_class, current_family, current_genus) {
    components <- strsplit(trimws(col), " ")[[1]]
    first_component <- components[1]
    
    result <- list(class = current_class, family = current_family, genus = current_genus, 
                   species = "", infraEpithet = "", scientificName = "")
    
    if (!grepl("[0-9]", first_component)) {
      # Class level
      result$class <- col
      result$family <- ""
      result$genus <- ""
    } else if (grepl("^[0-9]{2}$", first_component)) {
      # Family level
      result$family <- paste(components[-1], collapse = " ")
      result$genus <- ""
    } else if (grepl("^[0-9]{4}$", first_component)) {
      # Genus level
      genus <- components[2]
      # Replace "x" with "×" if present at the start of the genus
      if (startsWith(genus, "x")) {
        genus <- paste0("×", substr(genus, 2, nchar(genus)))
      }
      result$genus <- genus
    } else if (grepl("^[0-9]{6}[a-z]?$", first_component)) {
      # Species level
      full_name <- paste(components[-1], collapse = " ")
      
      # Replace abbreviated genus with full genus name
      if (grepl("^x?[A-Z]\\.", full_name)) {
        full_name <- sub("^x?[A-Z]\\.", result$genus, full_name)
      }
      
      result$species <- remove_infraEpithet(full_name)
      result$infraEpithet <- extract_infraEpithet(full_name)
      result$scientificName <- full_name
    }
    
    return(result)
  }
  
  for (i in 1:nrow(aba_selected)) {
    cat("\rSequencing name:", highcat(i), " / ", highcat(nrow(aba_selected)), fill = F)
    flush.console()
    
    row_result <- process_row(aba_selected$speciesCode[i], current_class, current_family, current_genus)
    aba_selected[i, c("class", "family", "genus", "species", "infraEpithet", "scientificName") := row_result]
    current_class <- row_result$class
    current_family <- row_result$family
    current_genus <- row_result$genus
  };catn()
  
  # Remove rows where all taxonomic fields are empty
  aba_selected <- aba_selected[!(scientificName == "")]
  
  aba_selected[, `:=`(
    class = clean_string(class, verbose),
    family = clean_string(family, verbose),
    genus = clean_string(genus, verbose),
    species = clean_string(species, verbose),
    infraEpithet = clean_string(infraEpithet, verbose),
    scientificName = clean_string(scientificName, verbose)
  )]
  
  # no need for speciesCode
  aba_selected[, speciesCode := NULL]
  
  vebcat("Formatting df.", veb = verbose)
  aba_formatted <- aba_selected
  
  # Remove all quotes in the names
  cols_to_clean <- c("class", "family", "genus", "species", "scientificName")
  
  ## Move the new columns to be the first four columns of the data table
  col_names <- names(aba_formatted)

  new_order <- c(col_names[1], col_names[(length(col_names) - 5):length(col_names)], col_names[2:(length(col_names) - 6)])
  
  aba_formatted <- aba_formatted[, ..new_order]
  ## Change to lowercase in class and family for cleaner look
  aba_formatted$class <- toupper(substr(aba_formatted$class, 1, 1)) %>% paste0(tolower(substr(aba_formatted$class, 2, nchar(aba_formatted$class))))
  aba_formatted$family <- toupper(substr(aba_formatted$family, 1, 1)) %>% paste0(tolower(substr(aba_formatted$family, 2, nchar(aba_formatted$family))))
  ## Replace the square symbol with hyphen in all columns
  aba_formatted[] <- lapply(aba_formatted, function(x) gsub("", "-", x))
  # Distinguish between present and absent data
  ## create a new dataset with species absent from the Arctic

  vebcat("Creating conditions.", veb = verbose)
  ## Create conditions for absence
  # Define columns to check
  cols_to_check <- names(aba_formatted)[which(names(aba_formatted) == "CW"):which(names(aba_formatted) == "YK")]
  aba_formatted[, borderline := borderline == 1]
  aba_formatted[, absent := all(.SD %in% c("-", "?", "**")), .SDcols = cols_to_check]

  vebcat("Making ABA present.", veb = verbose)
  aba_present <- aba_formatted[!(borderline | absent)]
  aba_present <- aba_present[, .(verbatimName, interimName = trimws(scientificName))]

  vebcat("Making ABA absent.", veb = verbose)
  aba_absent <- aba_formatted[borderline | absent]
  aba_absent <- aba_absent[, .(verbatimName, interimName = trimws(scientificName))]
  
  setnames(aba_present, "interimName", column)
  setnames(aba_absent, "interimName", column)
  
  lost <- nrow(aba_formatted) - (nrow(aba_present) + nrow(aba_absent))
  
  if (lost > 0) {
    vebcat("Error: formatted is not distributed equally in absent and present data tables.", color = "fatalError")
    stop("Something went wrong when wranling the ABA")
  }
  
  catn("Writing out files to:", colcat(dir, color = "output"))

  ?fwrite(aba_formatted, formatted_out, bom = TRUE)

  fwrite(aba_present, present_out, bom = TRUE)

  fwrite(aba_absent, absent_out, bom = TRUE)
  
  md_dt <- data.table(
    formatted = nrow(aba_formatted),
    present = nrow(aba_present),
    absent = nrow(aba_absent),
    lost = lost
  )

  mdwrite(
    config$files$post_seq_md,
    text = paste0("3;Wrangled ABA"),
    data = md_dt
  )

  return(list(
    present = aba_present,
    absent = aba_absent
  ))
}
