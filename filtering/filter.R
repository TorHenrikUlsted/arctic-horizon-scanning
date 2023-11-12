filtering = function() {
  message("------ Initiating filtering process ------ \n")
  
  ## Start filter script timer
  start_filter_time = Sys.time()
  
  source("./filtering/components/utils.R")
  source("./filtering/components/rgbif_similarity.R")
  source("./filtering/components/list_collector.R")
  
  components = collector()
  
  source("./filtering/components/synonym_check.R")
  
  # Ask if wanting to conduct synonym check
  synonym_check_input = readline(prompt = "Run synonym check (hours of waiting time)? 'n' loads pre-checked files. [y/n] ")
  
  # Source the synonym check if wanted
  if(tolower(synonym_check_input) == "y") {
    wfo_one_dfs = synonym_check(
      components$test_species,
      components$aba_present,
      components$aba_absent,
      components$ambio_present,
      components$ambio_absent,
      components$gbif_species,
      components$glonaf_species
      )
  } else {
    cat("Skipping synonym Check \n")
  }
  
  source("./filtering/components/list_merger.R")
  filtered_species = merger()
  
  a = ""
  while(a != "y" && a != "n") {
    a = readline(prompt = "Do you want to run a similarity check on the entire list (slow: ETC 30 min) [y], or do you want to check the 95% similar names (fast) [n]? ")
    
    if (a != "y" && a != "n") {
      cat("Invalid response. Please enter 'y' or 'n'.\n")
    }
  }
  ## Use the filtered species list with all the names
  if (a == "y") {
    scientific_names = unlist(filtered_species$scientificName)
    rgbif_simliarity_check(scientific_names)
  } 
  
  ## combine the two columns into one
  if (a == "n") {
    source("./filtering/components/jw_similarity.R")
    simCheck_filt_sp = jw_similarity_check(filtered_species)
    scientific_names = c(simCheck_filt_sp[ ,"name", drop = TRUE], simCheck_filt_sp[ ,"similarName", drop = TRUE])
    duplicate_species = scientific_names[!duplicated(scientific_names)]
    
    duplicate_species = data.frame(duplicate_species)
    names(duplicate_species)[names(duplicate_species) == 'duplicate_species'] = 'scientificName'
    filtered_species_jw = anti_join(filtered_species, duplicate_species, by = "scientificName") #Double check if this is the correct way to do it
    
    cat(yellow("Writing csv to: \n", "outputs/filtered_species_jw.csv \n"))
    write.csv(filtered_species_final, "outputs/filtered_species_jw.csv", row.names = F, fileEncoding = "UTF-8")
  }
  
  b = ""
  while(b != "y" && b != "n") {
    if (file.exists("./resources/gbif_occ_zip.zip") == T) {
      
    }
    a = readline(prompt = "Do you want to download GBIF occurences? ")
    
    if (b != "y" && b != "n") {
      cat("Invalid response. Please enter 'y' or 'n'.\n")
    }
  }
  
  
  # Stop the filter script timer
  end_filter_time = Sys.time()
  ### Calculate the time
  el_time_filter = format_elapsed_time(start_filter_time, end_filter_time)
  # Print message with elapsed time
  message("Script finished in ", el_time_filter)
}
