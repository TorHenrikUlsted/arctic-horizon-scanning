rgbif_simliarity_check <- function(scientific_names) {
  cat(blue("Initiating GBIF similarity check protocol. \n"))

  cat(cc$paleTurquoise("Using the name lookup method \n"))

  ## Only look in the vascular plants
  higher_taxon_key <- name_backbone(name = "Tracheophyta")$usageKey

  # Initialize variables to track progress
  ## pb = txtProgressBar(min = 0, max = length(scientific_names), style = 3, char="=")
  current_index <- 0

  speciesKeys <- lapply(scientific_names, function(x) {
    # Update the progress bar
    ## setTxtProgressBar(pb, which(scientific_names == x))
    # Update the current index
    current_index <<- current_index + 1

    # Display the current progress
    cat("\rProgress: ", current_index, " / ", length(scientific_names))

    # Flush the output buffer
    flush.console()

    # Measure the time it takes to process this scientific name
    result <- name_lookup(x, higherTaxonKey = higher_taxon_key)

    if ("speciesKey" %in% colnames(result$data)) {
      result$data$speciesKey
    } else {
      NA
    }
  })

  # Close the progress bar
  # close(pb)

  cat("\n Calculating and comparing species keys \n")
  ## Find the maximum number of speciesKey values
  max_speciesKeys <- max(sapply(speciesKeys, function(x) length(unique(x))))
  ## Create matrix
  mat <- matrix(nrow = length(speciesKeys), ncol = max_speciesKeys + 1)
  mat[, 1] <- scientific_names

  ## Error handling
  cat("Actual matrix dimensions: ", cc$paleTurquoise(dim(mat)), "\n")
  cat("Expected dimensions: ", cc$paleTurquoise(length(speciesKeys), max_speciesKeys + 1), "\n")

  dim_check <- dim(mat) - c(length(speciesKeys), max_speciesKeys + 1)
  if (all(dim_check == 0)) {
    cat("The dimensions of the matrix are ", green("correct"), "\n")
  } else {
    cat("The dimensions of the matrix are ", red("incorrect"), "\n")
  }

  # Add the species keys to the matrix
  for (i in seq_along(speciesKeys)) {
    unique_keys <- unique(speciesKeys[[i]])
    if (length(unique_keys) > 0) {
      mat[i, 2:(length(unique_keys) + 1)] <- unique_keys
    }
  }

  ## Convert matrix to data frame
  speciesKeys_unique <- as.data.frame(mat, stringAsFactors = FALSE)
  ## Convert the value columns to numeric
  speciesKeys_unique[, -1] <- lapply(speciesKeys_unique[, -1], as.numeric)
  ## Set the column names
  colnames(speciesKeys_unique) <- c("scientificName", "speciesKey")
  ## remove all other columns to only keep one value
  speciesKeys_unique <- speciesKeys_unique[, c(1, 2)]
  ## Get duplicated species
  duplicate_species <- speciesKeys_unique[duplicated(speciesKeys_unique$speciesKey), ]
  ## Filter duplicated species
  gbif_species_check <- speciesKeys_unique[!duplicated(speciesKeys_unique$speciesKey), ]

  cat(cc$lightSteelBlue("Filtering out", nrow(duplicate_species), "duplicated species based off of the rgbif speciesKeys \n"))

  cat(yellow("Wrinting duplicated species csv file to: \n", "outputs/similarity_check_outputs/similarityCheck_species_duplicated.csv \n"))

  create_dir_if("./outputs/filter/gbif-similarity-check")
  duplicate_species <- set_df_utf8(duplicate_species)
  fwrite(duplicate_species, "./outputs/filter/gbif-similarity-check/duplicated-species.csv", row.names = F, bom = T)
  
  check_savefile <- paste0("./outputs/filter/gbif-similarity-check/sim-checked-", region, ".csv")
  cat("Writing out similarity check csv to: \n", check_savefile, "\n")
  gbif_species_check <- set_df_utf8(gbif_species_check)
  fwrite(gbif_species_check, check_savefile, row.names = F, bom = T)

  cat(cc$lightGreen("GBIF similarity check protocol completed successfully. \n"))
  
  return(gbif_species_check)
}
