rgbif_simliarity_check = function(scientific_names) {
  message("------ Checking duplicates against GBIF database ------")
  
  cat("Using the name lookup method \n")
  
  ## Only look in the vascular plants
  higherTaxonKey = name_backbone(name = "Tracheophyta")$usageKey
  
  # Initialize variables to track progress
  ##pb = txtProgressBar(min = 0, max = length(scientific_names), style = 3, char="=")
  current_index = 0
  
  speciesKeys = lapply(scientific_names, function(x) {
    # Update the progress bar
    ##setTxtProgressBar(pb, which(scientific_names == x))
    # Update the current index
    current_index <<- current_index + 1
    
    # Display the current progress
    cat("\rProgress: ", current_index, " / ", length(scientific_names))
    
    # Flush the output buffer
    flush.console()
    
    # Measure the time it takes to process this scientific name
    result = name_lookup(x, higherTaxonKey = higherTaxonKey)
    
    if ("speciesKey" %in% colnames(result$data)) {
      result$data$speciesKey
    } else {
      NA
    }
  })
  
  # Close the progress bar
  #close(pb)
  
  cat("Calculating and comparing species keys \n")
  ## Find the maximum number of speciesKey values
  max_speciesKeys = max(sapply(speciesKeys, function(x) length(unique(x))))
  ## Create matrix
  mat = matrix(nrow = length(speciesKeys), ncol = max_speciesKeys + 1)
  mat[, 1] = scientific_names
  
  ## Error handling
  cat("Actual matrix dimensions: ", dim(mat), "\n")
  cat("Expected dimensions: ", length(speciesKeys), max_speciesKeys + 1, "\n")
  dim_check = dim(mat) - c(length(speciesKeys), max_speciesKeys + 1)
  if (all(dim_check == 0)) {
    cat("The dimensions of the matrix are correct \n")
  } else {
    cat("The dimensions of the matrix are incorrect \n")
  }
  
  # Add the species keys to the matrix
  for (i in seq_along(speciesKeys)) {
    unique_keys = unique(speciesKeys[[i]])
    if (length(unique_keys) > 0) {
      mat[i, 2:(length(unique_keys) + 1)] = unique_keys
    }
  }
  
  ## Convert matrix to data frame
  speciesKeys_unique = as.data.frame(mat, stringAsFactors = FALSE)
  ## Convert the value columns to numeric
  speciesKeys_unique[, -1] = lapply(speciesKeys_unique[, -1], as.numeric)
  ## Set the column names
  colnames(speciesKeys_unique) = c("scientificName", "speciesKey")
  ## remove all other columns to only keep one value
  speciesKeys_unique = speciesKeys_unique[, c(1, 2)]
  ## Get duplicated species
  duplicate_species = speciesKeys_unique[duplicated(speciesKeys_unique$speciesKey), ]
  ## Filter duplicated species
  filtered_duplicate_species = speciesKeys_unique[!duplicated(speciesKeys_unique$speciesKey), ]
  
  message("Filter out species")
  cat("Filtering out", nrow(duplicate_species), "Species from the filtered_species list: \n", paste(duplicate_species$scientificName, collapse = "\n"), "\n")
  
  write.csv(speciesKeys_unique, "outputs/similarity_check_outputs/similarityCheck_species_keys.csv", row.names = F, fileEncoding = "UTF-8")
  write.csv(filtered_duplicate_species, "outputs/similarity_check_outputs/filtered_duplicate_species.csv", row.names = F, fileEncoding = "UTF-8")
  
  return(
    filtered_duplicate_species
  )
}