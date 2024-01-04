synonym_check = function(test_species, aba_present, aba_absent, ambio_present, ambio_absent, gbif_species, glonaf_species) {
  
    simpleListNames = c("aba", "aba:present", "aba:absent", "ambio", "ambio:present", "ambio:absent", "gbif", "glonaf", "test")
    multiListNames = c("a", "n")
    inputString = ""
    wfoPromptMessage = 
      paste (
        "Which lists do you want to perform a synonym check on? 'n' loads pre-checked files. The possible commands are: \n list names: ", 
        paste(simpleListNames, collapse = ", "), 
        "\n multi-action(all, none): ", 
        paste(multiListNames, collapse = ", "), "\n"
      )
    inputString = checkInputCommands(inputString, simpleListNames, multiListNames, wfoPromptMessage)
    inputCommands = strsplit(inputString, ",")[[1]]
    inputCommands = trimws(inputCommands)
  
  test = function(test_species) {
    cat("Running the WFO synonym check for test species \n")
    cat("Number of species to analyse: ", yellow(nrow(test_species)), "\n")
    cat(yellow("Expected waiting time in hours: ", round((((nrow(test_species) * 8.45) / 60) / 60), digits = 2), "hours \n"))
    cat(yellow("Expected waiting time in minutes: ", round(((nrow(test_species) * 8.45) / 60), digits = 2), "minutes \n"))
    
    startTime = Sys.time()
    
    wfo_test_species = WFO.match(spec.data = test_species, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
    write.csv(wfo_test_species, "outputs/wfo_outputs/wfo_test_species.csv", row.names = F, fileEncoding = "UTF-8")
  
    wfo_one_test_species = WFO.one(WFO.result = wfo_test_species, priority = "Accepted", spec.name = "scientificName")
    write.csv(wfo_one_test_species, "outputs/wfo_one_outputs/wfo_one_test_species.csv", row.names = F, fileEncoding = "UTF-8")
    
    return(wfo_one_test_species)
    
    endTime = Sys.time()
    cat("WFO completed the match for test_species in ", format_elapsed_time(startTime, endTime), "\n")

  }
  
  aba_p = function(aba_present) {
    cat("Running the WFO synonym check for ABA present species \n")
    cat("Number of species to analyse: ", yellow(nrow(aba_present)), "\n")
    cat(yellow("Expected waiting time in hours: ", round((((nrow(aba_present) * 8.45) / 60) / 60), digits = 2), "hours \n"))
    cat(yellow("Expected waiting time in minutes: ", round(((nrow(aba_present) * 8.45) / 60), digits = 2), "minutes \n"))
      
    startTime = Sys.time()
      wfo_aba_arctic_present = WFO.match(spec.data = aba_present, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
      write.csv(wfo_aba_arctic_present, "outputs/wfo_outputs/wfo_aba_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      wfo_one_aba_arctic_present = WFO.one(WFO.result = wfo_aba_arctic_present, priority = "Accepted", spec.name = "scientificName")
      write.csv(wfo_one_aba_arctic_present, "outputs/wfo_one_outputs/wfo_one_aba_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(wfo_one_aba_arctic_present)
      
    endTime = Sys.time()
    cat("WFO completed the match for aba_present species in ", format_elapsed_time(startTime, endTime), "\n")
  }
    
  aba_a = function(aba_absent) {
    cat("Running the WFO synonym check for ABA absent species \n")
    cat("Number of species to analyse: ", yellow(nrow(aba_absent)), "\n")
    cat(yellow("Expected waiting time in hours: ", round((((nrow(aba_absent) * 8.45) / 60) / 60), digits = 2), "hours \n"))
    cat(yellow("Expected waiting time in minutes: ", round(((nrow(aba_absent) * 8.45) / 60), digits = 2), "minutes \n"))
    
      
    startTime = Sys.time()
      wfo_aba_arctic_absent = WFO.match(spec.data = aba_absent, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
      write.csv(wfo_aba_arctic_absent, "outputs/wfo_outputs/wfo_aba_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      wfo_one_aba_arctic_absent = WFO.one(WFO.result = wfo_aba_arctic_absent, priority = "Accepted", spec.name = "scientificName")
      write.csv(wfo_one_aba_arctic_absent, "outputs/wfo_one_outputs/wfo_one_aba_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(wfo_one_aba_arctic_absent)
      
    endTime = Sys.time()
    cat("WFO completed the match for aba_absent species in ", format_elapsed_time(startTime, endTime), "\n")
  }

    
  ambio_p = function(ambio_present) {
    cat("Running the WFO synonym check AMBIO present species \n")
    cat("Number of species to analyse: ", yellow(nrow(ambio_present)), "\n")
    cat(yellow("Expected waiting time in hours: ", round((((nrow(ambio_present) * 8.45) / 60) / 60), digits = 2), "hours \n"))
    cat(yellow("Expected waiting time in minutes: ", round(((nrow(ambio_present) * 8.45) / 60), digits = 2), "minutes \n"))
      
    startTime = Sys.time()

      wfo_ambio_arctic_present = WFO.match(spec.data = ambio_present, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
      write.csv(wfo_ambio_arctic_present, "outputs/wfo_outputs/wfo_ambio_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      wfo_one_ambio_arctic_present = WFO.one(WFO.result = wfo_ambio_arctic_present, priority = "Accepted", spec.name = "scientificName")
      write.csv(wfo_one_ambio_arctic_present, "outputs/wfo_one_outputs/wfo_one_ambio_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
    return(wfo_one_ambio_arctic_present)
      
    endTime = Sys.time()
    cat("WFO completed the match for ambio_present species in ", format_elapsed_time(startTime, endTime), "\n")
  }
    
  ambio_a = function(ambio_absent) {
    cat("Running the WFO synonym check for AMBIO absent species \n")
    cat("Number of species to analyse: ", yellow(nrow(ambio_absent)), "\n")
    cat(yellow("Expected waiting time in hours: ", round((((nrow(ambio_absent) * 8.45) / 60) / 60), digits = 2), "hours \n"))
    cat(yellow("Expected waiting time in minutes: ", round(((nrow(ambio_absent) * 8.45) / 60), digits = 2), "minutes \n"))
      
    startTime = Sys.time()
    wfo_ambio_arctic_absent = WFO.match(spec.data = ambio_absent, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
    write.csv(wfo_ambio_arctic_absent, "outputs/wfo_outputs/wfo_ambio_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
    
    wfo_one_ambio_arctic_absent = WFO.one(WFO.result = wfo_ambio_arctic_absent, priority = "Accepted", spec.name = "scientificName")
    write.csv(wfo_one_ambio_arctic_absent, "outputs/wfo_one_outputs/wfo_one_ambio_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
    return(wfo_one_ambio_arctic_absent)
    
    endTime = Sys.time()
      
    cat("WFO completed the match for ambio_absent species in ", format_elapsed_time(startTime, endTime), "\n")
  }

  
  gbif = function(gbif_species) {
    cat("Running the WFO synonym check for the GBIF list.\n")
    cat("Number of species to analyse: ", yellow(nrow(gbif_species)), "\n")
    cat(yellow("Expected waiting time in hours: ", round((((nrow(gbif_species) * 8.45) / 60) / 60), digits = 2), "hours \n"))
    cat(yellow("Expected waiting time in minutes: ", round(((nrow(gbif_species) * 8.45) / 60), digits = 2) , "minutes \n"))
    
    startTime = Sys.time()
      wfo_gbif_species = WFO.match(spec.data = gbif_species, spec.name = "gbif_species", WFO.file = WFO_file, verbose = T, counter = 1)
      write.csv(wfo_gbif_species, "outputs/wfo_outputs/wfo_gbif_species.csv", row.names = F, fileEncoding = "UTF-8")
      
      wfo_one_gbif_species = WFO.one(WFO.result = wfo_gbif_species, priority = "Accepted", spec.name = "scientificName")
      write.csv(wfo_one_gbif_species, "outputs/wfo_one_outputs/wfo_one_gbif_species.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(wfo_one_gbif_species)
    endTime = Sys.time()
    
    cat("WFO completed the match for the GBIF species in ", format_elapsed_time(startTime, endTime), "\n")
  }
  
  glonaf = function(glonaf_species) {
    cat("Running the WFO synonym check for the GloNAF species list.\n")
    cat("Number of species to analyse: ", yellow(nrow(glonaf_species)), "\n")
    cat(yellow("Expected waiting time in hours: ", round((((nrow(glonaf_species) * 8.45) / 60) / 60), digits = 2), "hours \n"))
    cat(yellow("Expected waiting time in minutes: ", round(((nrow(glonaf_species) * 8.45) / 60), digits = 2), "minutes \n"))
    
    startTime = Sys.time()
      wfo_glonaf_species = WFO.match(spec.data = glonaf_species, spec.name = "standardized_name", WFO.file = WFO_file, verbose = T, counter = 1)
      write.csv(wfo_glonaf_species, "outputs/wfo_outputs/wfo_glonaf_species.csv", row.names = F, fileEncoding = "UTF-8")
      
      wfo_one_glonaf_species = WFO.one(WFO.result = wfo_glonaf_species, priority = "Accepted", spec.name = "scientificName")
      write.csv(wfo_one_glonaf_species, "outputs/wfo_one_outputs/wfo_one_glonaf_species.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(wfo_one_glonaf_species)
    endTime = Sys.time()
    
    cat("WFO completed the match for the GloNAF species in ", format_elapsed_time(startTime, endTime), "\n")
  }
  
  
  if ("n" %in% inputCommands) {
    cat("None of the lists will be checked for synonyms. Moving on... \n")
  } else if ("a" %in% inputCommands) {
    aba_p(aba_present)
    aba_a(aba_absent)
    ambio_p(ambio_present)
    ambio_a(ambio_absent)
    gbif(gbif_species)
    glonaf(glonaf_species)

  } else {
    if ("test" %in% inputCommands) test(test_species)
    
    if ("aba" %in% inputCommands) {
      aba_p(aba_present)
      aba_a(aba_absent)

    }
    else if("aba:present" %in% inputCommands) aba_p(aba_present)
    else if("aba:absent" %in% inputCommands) aba_a(aba_absent)
    
    if ("ambio" %in% inputCommands) {
      ambio_p(ambio_present)
      ambio_a(ambio_absent)

    }
    else if("ambio:present" %in% inputCommands) ambio_p(ambio_present)
    else if("ambio:absent" %in% inputCommands) ambio_a(ambio_absent)
    
    if ("gbif" %in% inputCommands) gbif(gbif_species)
    if ("glonaf" %in% inputCommands) glonaf(gbif_species)
  }
  
  
}