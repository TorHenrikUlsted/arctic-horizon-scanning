rgbif_simliarity_check <- function(scientific_names, region.name, cores.max, taxon) {
  cat(blue("Initiating GBIF similarity check protocol. \n"))

  cat(cc$paleTurquoise("Using the name lookup method \n"))

  ## Only look in the vascular plants
  higher_taxon_key <- name_backbone(name = taxon)$usageKey

  cl <- makeCluster(cores.max)

  clusterEvalQ(cl, {
    library(rgbif)
  })
  
  csv_log_file <- "./outputs/filter/logs/sp-keys.csv"
  
  if (!file.exists(csv_log_file)) {
    file.create(csv_log_file)
  } else {
    file.remove(csv_log_file)
    file.create(csv_log_file)
  }
  
  log_file <- "./outputs/filter/logs/sp-similarity-progress.txt"
  
  if (!file.exists(log_file)) {
    file.create(log_file)
  } else {
    file.remove(log_file)
    file.create(log_file)
  }
  
  try(log_file <- file(log_file, open = "at"))
  sink(log_file, type = "output")
  

  gbif_sim_timer <- start_timer("gbif_names_timer")

  sp_keys <- clusterApplyLB(cl, seq_along(scientific_names), function(i) {
    tryCatch(
      {

        cat("collecting species\n",
            sprintf("%6s | %8s | %7s\n", "key", "Max keys", "Percent"),
            sprintf("\r%6d | %8d | %7.2f %%", i, length(scientific_names), round(i / length(scientific_names) * 100, 2)))
        
        result <- name_lookup(i, higherTaxonKey = higher_taxon_key)

        speciesKey <- if ("speciesKey" %in% colnames(result$data)) {
          result$data$speciesKey[1]
        } else {
          NA
        }
        
        speciesKey_df <- data.frame(
          iteration = i,
          speciesKey = speciesKey
        )
        
        # Write the result and the iteration number to a file
        fwrite(speciesKey_df, file = "./outputs/filter/logs/sp-keys.csv", append = TRUE)
        
        invisible(gc())
        
        return(speciesKey)
        
      },
      error = function(e) {
        e <- conditionMessage(e)
        print(e)
      }
    )
  })
  
  sink(type = "output")
  close(log_file)


  cat("Finishing up \n")

  stopCluster(cl)

  end_timer(gbif_sim_timer)

  cat(cc$lightSteelBlue(length(sp_keys)), "GBIF species keys found. \n")

  cat("\n Calculating and comparing species keys \n")
  ## Find the maximum number of speciesKey values
  max_speciesKeys <- max(sapply(sp_keys, function(x) length(unique(x))))
  ## Create matrix
  mat <- matrix(nrow = length(sp_keys), ncol = max_speciesKeys + 1)
  mat[, 1] <- scientific_names

  ## Error handling
  cat("Actual matrix dimensions: ", cc$paleTurquoise(dim(mat)), "\n")
  cat("Expected dimensions: ", cc$paleTurquoise(length(sp_keys), max_speciesKeys + 1), "\n")

  dim_check <- dim(mat) - c(length(sp_keys), max_speciesKeys + 1)
  if (all(dim_check == 0)) {
    cat("The dimensions of the matrix are ", green("correct"), "\n")
  } else {
    cat("The dimensions of the matrix are ", red("incorrect"), "\n")
  }

  # Add the species keys to the matrix
  for (i in seq_along(sp_keys)) {
    unique_keys <- unique(sp_keys[[i]])
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

  check_dups <- paste0("./outputs/filter/gbif-acquisition/gbif-duplicated-", region.name, ".csv")
  cat(yellow("Writing duplicated species csv file to: \n", check_dups, "\n"))
  create_dir_if("./outputs/filter/gbif-acquisition")
  duplicate_species <- set_df_utf8(duplicate_species)
  fwrite(duplicate_species, check_dups, row.names = F, bom = T)

  check_savefile <- paste0("./outputs/filter/gbif-acquisition/gbif-sim-checked-", region.name, ".csv")
  cat("Writing out similarity check csv to: \n", check_savefile, "\n")
  gbif_species_check <- set_df_utf8(gbif_species_check)
  fwrite(gbif_species_check, check_savefile, row.names = F, bom = T)

  cat(cc$lightGreen("GBIF similarity check protocol completed successfully. \n"))

  return(gbif_species_check)
}
