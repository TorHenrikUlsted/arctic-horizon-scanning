synonym_check =  function() {
  
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
  
  
  test = function() {
      message("Running WFO Synonym Check [TEST]... \n")
      Sys.sleep(2)
  }
  
  aba_p = function() {
    cat("Running the WFO synonym check for ABA present species \n")
    cat("Number of species to analyse: ", nrow(aba_present), "\n")
      
    startTime = Sys.time()
      wfo_aba_arctic_present = WFO.match(spec.data = aba_present, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
      wfo_one_aba_arctic_present = WFO.one(WFO.result = wfo_aba_arctic_present, priority = "Accepted", spec.name = "Species")
    endTime = Sys.time()
      
    cat("WFO completed the match for aba_present species in ", format_elapsed_time(startTime, endTime), "\n")
    cat("Creating CSV files of the outputs \n")
      write.csv(wfo_aba_arctic_present, "outputs/wfo_outputs/wfo_aba_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      write.csv(wfo_one_aba_arctic_present, "outputs/wfo_one_outputs/wfo_one_aba_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
  }
    
  aba_a = function() {
    cat("Running the WFO synonym check for ABA absent species \n")
    cat("Number of species to analyse: ",nrow(aba_absent), "\n")
      
    startTime = Sys.time()
      wfo_aba_arctic_absent = WFO.match(spec.data = aba_absent, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
      wfo_one_aba_arctic_absent = WFO.one(WFO.result = wfo_aba_arctic_absent, priority = "Accepted", spec.name = "Species")
    endTime = Sys.time()
      
    cat("WFO completed the match for aba_absent species in ", format_elapsed_time(startTime, endTime), "\n")
    cat("Creating CSV file of the output \n")
      write.csv(wfo_aba_arctic_absent, "outputs/wfo_outputs/wfo_aba_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      write.csv(wfo_one_aba_arctic_absent, "outputs/wfo_one_outputs/wfo_one_aba_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
  }

  
    
  ambio_p = function() {
    cat("Running the WFO synonym check AMBIO present species \n")
    cat("Number of species to analyse: ",nrow(ambio_present), "\n")
      
    startTime = Sys.time()
      wfo_ambio_arctic_present = WFO.match(spec.data = ambio_present, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
      wfo_one_ambio_arctic_present = WFO.one(WFO.result = wfo_ambio_arctic_present, priority = "Accepted", spec.name = "Species")
    endTime = Sys.time()
      
    cat("WFO completed the match for ambio_present species in ", format_elapsed_time(startTime, endTime), "\n")
    cat("Creating CSV file of the output \n")
      write.csv(wfo_ambio_arctic_present, "outputs/wfo_outputs/wfo_ambio_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      write.csv(wfo_one_ambio_arctic_present, "outputs/wfo_one_outputs/wfo_one_ambio_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
  }
    
  ambio_a = function() {
    cat("Running the WFO synonym check for AMBIO absent species \n")
    cat("Number of species to analyse: ",nrow(ambio_absent), "\n")
      
    startTime = Sys.time()
    wfo_ambio_arctic_absent = WFO.match(spec.data = ambio_absent, spec.name = "Species_SubSpecies", WFO.file = WFO_file, verbose = T, counter = 1)
    wfo_one_ambio_arctic_absent = WFO.one(WFO.result = wfo_ambio_arctic_absent, priority = "Accepted", spec.name = "Species")
    endTime = Sys.time()
      
    cat("WFO completed the match for ambio_absent species in ", format_elapsed_time(startTime, endTime), "\n")
      
    cat("Creating CSV file of the output \n")
      write.csv(wfo_ambio_arctic_absent, "outputs/wfo_outputs/wfo_ambio_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      write.csv(wfo_one_ambio_arctic_absent, "outputs/wfo_one_outputs/wfo_one_ambio_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
  }

  
  gbif = function() {
    cat("Running the WFO synonym check for the GBIF list.\n")
    cat("Number of species to analyse: ",nrow(gbif_species), "\n")
    
    startTime = Sys.time()
      wfo_gbif_species = WFO.match(spec.data = gbif_species, spec.name = "gbif_species", WFO.file = WFO_file, verbose = T, counter = 1)
      wfo_one_gbif_species = WFO.one(WFO.result = wfo_gbif_species, priority = "Accepted", spec.name = "Species")
    endTime = Sys.time()
    
    cat("WFO completed the match for the GBIF species in ", format_elapsed_time(startTime, endTime), "\n")
    
    cat("Creating CSV file of the output \n")
      write.csv(wfo_gbif_species, "outputs/wfo_outputs/wfo_gbif_species.csv", row.names = F, fileEncoding = "UTF-8")
      write.csv(wfo_one_gbif_species, "outputs/wfo_one_outputs/wfo_one_gbif_species.csv", row.names = F, fileEncoding = "UTF-8")
  }
  
  glonaf = function() {
    cat("Running the WFO synonym check for the GloNAF species list.\n")
    cat("Number of species to analyse: ",nrow(glonaf_species), "\n")
    
    startTime = Sys.time()
      wfo_glonaf_species = WFO.match(spec.data = glonaf_species, spec.name = "standardized_name", WFO.file = WFO_file, verbose = T, counter = 1)
      wfo_one_glonaf_species = WFO.one(WFo.result = wfo_glonaf_species, priority = "Accepted", spec.name = "Species")
    endTime = Sys.time()
    
    cat("WFO completed the match for the GloNAF species in ", format_elapsed_time(startTime, endTime), "\n")

    cat("Creating CSV file of the output \n")
      write.csv(wfo_glonaf_species, "outputs/wfo_outputs/wfo_glonaf_species.csv", row.names = F, fileEncoding = "UTF-8")
      write.csv(wfo_one_glonaf_species, "outputs/wfo_one_outputs/wfo_one_glonaf_species.csv", row.names = F, fileEncoding = "UTF-8")
  }
  
  
  if ("n" %in% inputCommands) {
    cat("None of the lists will be checked for synonyms. Moving on... \n")
  } else if ("a" %in% inputCommands) {
    aba_p()
    aba_a()
    ambio_p()
    ambio_a()
    gbif()
    glonaf()
    
  } else {
    
    if ("test" %in% inputCommands) test()
    
    if ("aba" %in% inputCommands) {
      aba_p()
      aba_a()
    }
    else if("aba:present" %in% inputCommands) aba_p()
    else if("aba:absent" %in% inputCommands) aba_a()
    
    if ("ambio" %in% inputCommands) {
      ambio_p()
      ambio_a()
    }
    else if("ambio:present" %in% inputCommands) ambio_p()
    else if("ambio:absent" %in% inputCommands) ambio_a()
    
    if ("gbif" %in% inputCommands) gbif()
    if ("glonaf" %in% inputCommands) glonaf()
  }
  
}
