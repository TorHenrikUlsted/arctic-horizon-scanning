# Filter all species lists

initiate = function() {
  ## Start filter script timer
  start_filterScript_time <<- Sys.time()
  ## Source scripts ETC 10 min without the synonym check
  startTime <<- Sys.time()
  ## Source utility script
  source("components/utils.R")
  ## create the possible commands
  simpleListNames <<- c("aba", "ambio", "gbif", "glonaf", "test")
  multiListNames <<- c("all", "none")
  ## create an empty string of the input
  inputString = ""
  ## specify the prompt message
  wfoPromptMessage = paste("Which lists do you want to run a synonym check on? The possible commands are: \n list names: ", paste(simpleListNames, collapse = ", "), "\n multi-action: ", paste(multiListNames, collapse = ", "), "\n")
  ## Run the command line input check
  inputString = checkInputCommands(inputString, simpleListNames, multiListNames, wfoPromptMessage)
  ## Split input into individual commands
  inputCommands <<- strsplit(inputString, ",")[[1]]
  inputCommands <<- trimws(inputCommands)
  
  message("Sourcing components. \n")
  sourcing()
} 

sourcing = function() {
  ## Define components and their associated file paths
  components = list(
    aba = list(
      source = "components/aba_wrangling.R",
      present = "outputs/wfo_one_aba_arctic_present.csv",
      absent = "outputs/wfo_one_aba_arctic_absent.csv"
    ),
    ambio = list(
      source = "components/ambio_wrangling.R",
      present = "outputs/wfo_one_ambio_arctic_present.csv",
      absent = "outputs/wfo_one_ambio_arctic_absent.csv"
    ),
    gbif = list(
      source = c("components/regions_data_import.R", "components/gbif_species_retrieval.R"),
      species = "outputs/wfo_one_gbif_species.csv"
    ),
    glonaf = list(
      source = "components/glonaf_wrangling.R",
      species = "outputs/wfo_one_glonaf_species.csv"
    )
  )
  
  ## Handle "all" and "none" commands
  if ("all" %in% inputCommands) {
    inputCommands = names(components)
  } else if ("none" %in% inputCommands) {
    inputCommands = character(0)
  }
  
  ## Iterate over components
  for (component in names(components)) {
    if (component %in% inputCommands) {
      ## Source component files
      for (file in components[[component]]$source) {
        source(file)
      }
    } else {
      ## Load data from files
      for (field in setdiff(names(components[[component]]), "source")) {
        file <- components[[component]][[field]]
        tryCatch({
          var_name <- paste0("wfo_", component, "_", field)
          assign(var_name, read.csv(file), envir = .GlobalEnv)
        }, warning = function(w) {
          print(paste("Warning:", w))
        }, error = function(e) {
          print(paste("Error:", e))
        })
      }
    }
  }
  
  # Ask if wanting to conduct synonym check
  synonym_check_input = readline(prompt = "Run synonym check (hours of waiting time)? [y/n] ")
  
  # Source the synonym check if wanted
  if(tolower(synonym_check_input) == "y") {
    source("components/synonym_check.R")
  } else {
    cat("Skipping synonym Check \n")
  }
}

message("Initiating filtering process \n")
initiate()

# Stop sourcing timer
endTime = Sys.time()
### Calculate the time
formatted_elapsed_time = format_elapsed_time(startTime, endTime)
## Print message with elapsed time
message("Sourcing script finished in ", formatted_elapsed_time)

remove_duplicates = function() {
  # choose Scientifically accepted names
  message("---------- Selecting and removing duplicates from the ABA, AMBIO, GBIF, and GloNAF lists ----------")
  
  ## Define components and their associated file paths
  components = list(
    aba = list(
      present = "wfo_one_aba_present",
      absent = "wfo_one_aba_absent"
    ),
    ambio = list(
      present = "wfo_one_ambio_present",
      absent = "wfo_one_ambio_absent"
    ),
    gbif = list(species = "wfo_one_gbif_species"),
    glonaf = list(species = "wfo_one_glonaf_species")
  )
  
  ## Iterate over components
  for (component in names(components)) {
    if (component != "gbif" && component != "glonaf") {
      present_data = get(components[[component]]$present)
      absent_data = get(components[[component]]$absent)
      
      ## Select the scientific names and remove NA
      present_data = present_data %>% 
        select(scientificName) %>% 
        filter(!is.na(scientificName))
      cat(paste0("All ", component, "_arctic_present NAs removed: ", any(!is.na(present_data)), "\n"))
      
      ## Remove duplicates
      present_data = distinct(present_data)
      cat(paste0("All ", component, "_arctic_present are distinct: ", any(!duplicated(present_data)), "\n"))
      
      ## Select the scientific names and remove NA
      absent_data = absent_data %>% 
        select(scientificName) %>% 
        filter(!is.na(scientificName))
      cat(paste0("All ", component, "_arctic_absent NAs removed: ", any(!is.na(absent_data)), "\n"))
      
      ## Remove duplicates
      absent_data = distinct(absent_data)
      cat(paste0("All ", component, "_arctic_absent are distinct: ", any(!duplicated(absent_data)), "\n"))
      
      ## Assign data to global environment
      assign(paste0(component, "_arctic_present"), present_data, envir = .GlobalEnv)
      assign(paste0(component, "_arctic_absent"), absent_data, envir = .GlobalEnv)
    } else {
      species_data = get(components[[component]]$species)
      species_data = species_data %>% 
        select(scientificName) %>% 
        filter(!is.na(scientificName))
      
      ## Remove duplicates
      species_data = distinct(species_data)
      
      ## Assign data to global environment
      assign(paste0(component, "_species"), species_data, envir = .GlobalEnv)
    }
  }
  
  
  # Combine AMBIO and ABA lists
  message("------ Merging ABA and AMBIO lists ------")
  
  ## merge present lists
  arctic_present = aba_arctic_present %>% 
    full_join(ambio_arctic_present, by = "scientificName") %>% 
    distinct(scientificName, .keep_all = TRUE)
  
  ## Run NA and distinct check
  cat("All arctic_present NAs removed:", any(!is.na(arctic_present)), "\n")
  cat("All arctic_present values are distinct: ", any(!duplicated(arctic_present)), "\n")
  
  ## merge absent lists
  arctic_absent = aba_arctic_absent %>% 
    full_join(ambio_arctic_absent, by = "scientificName") %>% 
    distinct(scientificName, .keep_all = TRUE)
  
  ## Run NA and distinct check
  cat("All arctic_absent NAs removed:", any(!is.na(arctic_absent)), "\n")
  cat("All arctic_absent values are distinct: ", any(!duplicated(arctic_absent)), "\n")
  
  # Merge GBIF list with ABA and AMBIO list
  message("------ Merging the combined ABA and AMBIO list with the GBIF list ------")
  
  ## Remove rows from filtered_species data frame that have matching values to the arctic_absent species. This is to avoid duplicates.
  gbif_arctic_present = merge(gbif_species, 
                              arctic_present, 
                              by = "scientificName"
  )
  ## Remove rows from gbif_species data frame that have matching values to the arctic_present species
  gbif_arctic_absent = anti_join(gbif_species, 
                                 arctic_present, 
                                 by = "scientificName"
  )
  
  ## Remove rows that are similar to the arctic_absent to avoid duplicates
  gbif_arctic_absent = anti_join(gbif_arctic_absent, 
                                 arctic_absent, 
                                 by = "scientificName"
  )
  
  ## Check how many are similar
  cat("Number of rows in GBIF that are similar to the arctic_absent list: ", nrow(anti_join(gbif_species, arctic_present, by = "scientificName")) - nrow(anti_join(gbif_arctic_absent, arctic_absent, by = "scientificName")), "\n")
  
  ## Now add all the absent species to the list. 
  gbif_arctic_absent = bind_rows(gbif_arctic_absent, arctic_absent)
  
  cat("gbif_arctic_present and gbif_arctic_absent lists have been made\n")
  
  # Merge GloNAF list with ABA, AMBIO, and GBIF list
  message("------ Merging the combined ABA, AMBIO, and GBIF list with the GloNAF list ------")
  
  ## Remove rows from glonaf_species data frame that have matching values to the arctic_present species
  glonaf_arctic_present = merge(glonaf_species, 
                                arctic_present, 
                                by = "scientificName"
  )
  
  ## Remove rows from glonaf_species data frame that have matching values to the arctic_present species
  ## This one will merge with arctic_absent species later
  glonaf_arctic_absent = anti_join(glonaf_species, 
                                   arctic_present, 
                                   by = "scientificName"
  )
  
  cat("All glonaf_arctic_absent NAs removed:", any(!is.na(glonaf_arctic_absent)))
  ## Remove duplicates
  glonaf_arctic_absent = distinct(glonaf_arctic_absent)
  cat("All glonaf_arctic_absent are distinct:", any(!duplicated(glonaf_arctic_absent)))
  
  # Merge the filtered list with the GloNAF list
  message("------ Merging the gbif_arctic_absent list with the glonaf_arctic_absent list ------")
  ## Remove the glonaf species from the filtered_list. This is to avoid duplicates when merging
  filtered_species = anti_join(glonaf_arctic_absent, 
                               gbif_arctic_absent,
                               by = "scientificName"
  )
  
  ## Check how many are similar
  cat("Number of rows in gbif_arctic_absent that are similar to the glonaf_arctic_absent list: ", nrow(merge(glonaf_arctic_absent, gbif_arctic_absent)), "\n")
  
  ## Merge filtered_species from glonaf with the gbif_species
  filtered_species = bind_rows(gbif_arctic_absent, filtered_species)
  
  ## Check for NAs and duplications
  cat("All filtered_species NAs removed:", any(!is.na(filtered_species)), "\n")
  cat("All filtered_species are distinct:", any(!duplicated(filtered_species)), "\n")
  
  ## Replace odd multiplication sign with x in order to better display the final list
  filtered_species$scientificName = lapply(filtered_species$scientificName, function(x) gsub("(\\w) × (\\w)", "\\1 x \\2", x))
  filtered_species$scientificName = unlist(filtered_species$scientificName)
}


message("Removing duplicates")
remove_duplicates()

check_similarities = function() {
  ## Run a similarity check to find possible missed duplications
  simCheck_filt_sp = similarity_check(filtered_species, "scientificName", "scientificName", "jw", 0.05)
  simCheck_filt_sp = data.frame(simCheck_filt_sp)
  
  ## Filter out different variances
  simCheck_filt_sp = filter_rows_after_split_text(simCheck_filt_sp, "name", "similarName", "var.")
  ## Filter out hybrids
  simCheck_filt_sp = filter_rows_around_split_text(simCheck_filt_sp, "name", "similarName", "x")
  ## Filter out different subspecies
  simCheck_filt_sp = filter_rows_after_split_text(simCheck_filt_sp, "name", "similarName", "subsp.")
  ## Filter out only Genus names
  
  
  ## Further use rgbif to check the similar names
  rgif_check = function() {
    ## Use rgbif to check the similar names
    ## combine the two columns into one
    scientific_names = unlist(simCheck_filt_sp[ ,c("name", "similarName")] )
    
    ## Only look in the vascular plants
    higherTaxonKey = name_backbone(name = "Tracheophyta")$usageKey
    
    speciesKeys = lapply(scientific_names, function(x) {
      name_lookup(x, higherTaxonKey = higherTaxonKey)$data$speciesKey
    })
    
    ## Find the maximum number of speciesKey values
    max_speciesKeys = max(sapply(speciesKeys, function(x) length(unique(x))))
    ## Create matrix
    mat = matrix(nrow = length(speciesKeys), ncol = max_speciesKeys + 1)
    mat[, 1] = scientific_names
    for (i in seq_along(speciesKeys)) {
      mat[i, 2:(length(unique(speciesKeys[[i]])) + 1)] = unique(speciesKeys[[i]])
    }
    ## Convert matrix to data frame
    speciesKeys_unique = as.data.frame(mat, stringAsFactors = FALSE)
    ## Convert the value columns to numeric
    speciesKeys_unique[, -1] = lapply(speciesKeys_unique[, -1], as.numeric)
    ## Set the column names
    colnames(speciesKeys_unique) = c("scientificName", "speciesKey")
    ## remove all other columns to only keep one value
    speciesKeys_unique = speciesKeys_unique[, c(1, 2)]
    ## Write it out
    write.csv(speciesKeys_unique, "outputs/similarityCheck_filtered_species.csv", fileEncoding = "UTF-8")
    ## Get duplicated species
    duplicate_species = speciesKeys_unique[duplicated(speciesKeys_unique$speciesKey), ]
    ## Filter duplicated species
    filtered_duplicate_species = speciesKeys_unique[!duplicated(speciesKeys_unique$speciesKey), ]
    
    cat("Filtering out", nrow(duplicate_species), "Species from the filtered_species list: \n", paste(duplicate_species$scientificName, collapse = "\n"), "\n")
    
    ## Make the final list
    filtered_species_final = anti_join(filtered_species, duplicate_species, by = "scientificName")
    
  }
  message("Checking duplicates against GBIF database")
  rgif_check()
}
message("Running similarity check")
check_similarities()

write_lists = function() {
  ## Write it out
  write.csv(speciesKeys_unique, "outputs/similarityCheck_filtered_species.csv", fileEncoding = "UTF-8")

  ## Write out the final CSV of all species outside the Arctic
  write.csv(filtered_species_final, "outputs/filtered_species_final.csv", fileEncoding = "UTF-8")
}

message("writing out lists")
write_lists()

# Stop the filter script timer
end_filterScript_time = Sys.time()
### Calculate the time
el_time_filter = format_elapsed_time(start_filterScript_time, end_filterScript_time)
# Print message with elapsed time
message("Script finished in ", el_time_filter, "\n")