collector = function() {
  ## Source scripts ETC 10 min without the synonym check
  startTime = Sys.time()
  
  
  ## Define components and their associated file paths
  components = list(
    aba = list(
      source = "components/aba_wrangling.R",
      present = "outputs/wrangling_outputs/aba/aba_arctic_present.csv",
      absent = "outputs/wrangling_outputs/aba/aba_arctic_absent.csv"
    ),
    ambio = list(
      source = "components/ambio_wrangling.R",
      present = "outputs/wrangling_outputs/ambio/ambio_arctic_present.csv",
      absent = "outputs/wrangling_outputs/ambio/ambio_arctic_absent.csv"
    ),
    gbif = list(
      source = c("components/regions_data_import.R", "components/gbif_species_retrieval.R"),
      species = "outputs/gbif_retrieval_process/gbif_species_check.csv"
      #present = "outputs/gbif_retrieval_process/gbif_cavm_species.csv",
      #absent = "outputs/gbif_retrieval_process/gbif_boreal_species.csv"
    ),
    glonaf = list(
      source = "components/glonaf_wrangling.R",
      species = "outputs/wrangling_outputs/glonaf/glonaf_species_names.csv"
    ),
    test = list(
      source = "components/test.R",
      species = "outputs/wrangling_outputs/test/test_species.csv"
    ) 
  )
  
  
  missing_source_files = character(0)  # Initialize an empty vector to store missing source file names
  missing_check_files = character(0)   # Initialize an empty vector to store components that need to be checked
  
  # Iterate through components to check for missing files
  for (component_name in names(components)) {
    component = components[[component_name]]
    
    # Check source file
    if (!all(file.exists(component$source))) {
      missing_source_files = c(missing_source_files, component$source)
    }
    
    # Check other files (species, present, absent)
    other_files = c(component$species, component$present, component$absent)
    missing_other_files_component = other_files[!file.exists(other_files)]
    
    if (length(missing_other_files_component) > 0) {
      missing_check_files = c(missing_check_files, component_name)
    }
  }
  
  if (length(missing_source_files) > 0) {
    cat(red("The following files are missing and needed in order to continue:\n"))
    for (file_path in missing_source_files) {
      cat(cyan("- ", file_path, "\n"))
    }
  }
  
  if (length(missing_check_files) > 0) {
    cat(yellow("The following components are not wrangled and need to be run:\n"))
    for (component_name in missing_check_files) {
      cat(cyan("- ", component_name, "\n"))
    }
  }
  
  if (length(missing_source_files) == 0 && length(missing_check_files) == 0) {
    cat(green("All files needed to run the script already exist. No need to run any extra checks. \n", "Run lists you have edited \n"))
  }
  
  
  # Create command line
  simpleListNames = c("aba", "ambio", "gbif", "glonaf", "test")
  multiListNames = c("a", "n")
  inputString = ""
  wranglePromptMessage = paste("Which lists do you want to wrangle? load pre-wrangled files with 'n'. The possible commands are: \n list names: ", 
                               paste(simpleListNames, collapse = ", "), 
                               "\n multi-action(all, none): ", 
                               paste(multiListNames, collapse = ", "), 
                               "\n"
  )
  
  inputString = checkInputCommands(inputString, simpleListNames, multiListNames, wranglePromptMessage)
  inputCommands = strsplit(inputString, ",")[[1]]
  inputCommands = trimws(inputCommands)
  
  
  ## Handle "all" and "none" commands
  if ("a" %in% inputCommands) {
    inputCommands = names(components)
  } else if ("n" %in% inputCommands) {
    inputCommands = character(0)
  }
  
  # Create an empty list to store the data frames
  data_frames = list()
  
  if ("gbif" %in% inputCommands) {
    for (component in components[["gbif"]]$source) {
      source(component)
    }
    species_occ = sp_retrieve()
    data_frames[["species_occ"]] = species_occ
  }
  
  for (component in names(components)) {
    if (component %in% inputCommands) {
      ## Source component files
      for (file in components[[component]]$source) {
        source(file)
      }
    } else {
      for (field in setdiff(names(components[[component]]), "source")) {
        file = components[[component]][[field]]
        tryCatch({
          var_name = paste0(component, "_", field)
          data_frame = read.csv(file)
          # Store the data frame in the list
          data_frames[[var_name]] = data_frame
        }, warning = function(w) {
          print(paste("Warning:", w))
        }, error = function(e) {
          print(paste("Error:", e))
        })
      }
    }
  }
  

  
  # Return the list of data frames
  return(data_frames)
  
  # Stop sourcing timer
  endTime = Sys.time()
  ### Calculate the time
  formatted_elapsed_time = format_elapsed_time(startTime, endTime)
  ## Print message with elapsed time
  message("Sourcing script finished in ", formatted_elapsed_time)
}
