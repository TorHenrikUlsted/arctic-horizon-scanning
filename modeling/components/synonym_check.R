# Check for synonyms using the World Flora Online package
## create the possible commands
local_validCommands = c("aba", "ambio", "gbif", "glonaf")
global_validCommands = c("all", "none")
## create an empty string of the input
inputString <<- ""
## specify the prompt message
promptMessage = paste("Which lists do you want to run a synonym check on? The possible commands are: \n Single lists: ", paste(local_validCommands, collapse = ", "), "\n Multilists: ", paste(global_validCommands, collapse = ", "), "\n")
## Run the command line input check
inputString = checkInputCommands(inputString, local_validCommands, global_validCommands, promptMessage)

## Split input into individual commands
inputCommands = strsplit(inputString, ",")[[1]]


##check typed in commands and execute
if (any(inputCommands %in% global_validCommands || inputCommands %in% local_validCommands)) { # main if statement
       
    if ("none" %in% inputCommands) {
    cat("None of the lists will be checked for synonyms. Moving on... \n")
    
  } else { #start of all list if statements
    
    # -------------------- ABA -------------------- #

    if ("aba" %in% inputCommands || "all" %in% inputCommands) {
      cat("Running the WFO synonym check for the ABA list.\n")
      # Conduct WFO synonym matching
      ## Synonym check for the Arctic present species
      cat("Starting synonym check for ABA present species \n")
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
      
      ## Synonym check for Arctic absent species
      cat("Starting synonym check for ABA absent species \n")
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
      
    } else {
      ## Do not run the synonym check for this list
      cat("Skipping the WFO synonym check for the ABA list. \n")
    }
    
    # -------------------- AMBIO -------------------- #
    
    if ("ambio" %in% inputCommands || "all" %in% inputCommands) {
      cat("Running the WFO synonym check for the AMBIO list.\n")

        ## Synonym check for the Arctic present species
      cat("Starting synonym check for AMBIO present species \n")
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
        
        
        ## Synonym check for Arctic absent species
        cat("Starting synonym check for AMBIO absent species \n")
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
        
      
    } else {
      cat("Skipping the WFO synonym check for the AMBIO list. \n")
    }
    
    # -------------------- GBIF -------------------- #
    
    if ("gbif" %in% inputCommands || "all" %in% inputCommands) {
      cat("Running the WFO synonym check for the GBIF list.\n")
      
      ## Synonym check for GBIF
      ### Take the time
      start_time = Sys.time()
      wfo_gbif_species = WFO.match(spec.data = gbif_species, spec.name = "gbif_species", WFO.file = WFO_file, verbose = T, counter = 500)
      end_time = Sys.time()
      ### Calculate the time
      formatted_elapsed_time = format_elapsed_time(start_time, end_time)
      ## Print message with elapsed time
      cat("WFO completed the match for the GBIF species in ", formatted_elapsed_time, "\n")
      ## write it into a CSV file
      cat("Creating CSV file of the output \n")
      write.csv(wfo_gbif_species, "outputs/wfo_gbif_species.csv")
      
    } else {
      cat("Skipping the WFO synonym check for the GBIF list. \n")
    } 
    
    # -------------------- GloNAF -------------------- #
    
    if ("glonaf" %in% inputCommands || "all" %in% inputCommands) {
      cat("Running the WFO synonym check for the GloNAF species list.\n")
      
      ## Synonym check for GloNAF
      ### Take the time
      start_time = Sys.time()
      wfo_glonaf_species = WFO.match(spec.data = glonaf_species, spec.name = "standardized_name", WFO.file = WFO_file, verbose = T, counter = 500)
      end_time = Sys.time()
      ### Calculate the time
      formatted_elapsed_time = format_elapsed_time(start_time, end_time)
      ## Print message with elapsed time
      cat("WFO completed the match for the GloNAF species in ", formatted_elapsed_time, "\n")
      ## write it into a CSV file
      cat("Creating CSV file of the output \n")
      write.csv(wfo_glonaf_species, "outputs/wfo_glonaf_species.csv")
      
    } else {
      cat("Skipping the WFO synonym check for the GloNAF list. \n")
    }
    cat("Finished the synonym checks")
  } # end of all list if statements
} # end of main if statement