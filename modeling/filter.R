filtering = function() {
  message("------ Initiating filtering process ------ \n")
  
  ## Start filter script timer
  start_filter_time = Sys.time()
  
  ## Get and Set working directory
  tryCatch({
    setwd("modeling")
  }, error = function(e) {
    cat("Working directory: ", getwd(), "\n")
  })
  
  getwd()
  source("components/utils.R")
  source("components/list_collector.R")
  source("components/synonym_check.R")
  source("components/list_merger.R")
  
  collector()
  
  # Ask if wanting to conduct synonym check
  synonym_check_input = readline(prompt = "Run synonym check (hours of waiting time)? 'n' loads pre-checked files. [y/n] ")
  
  # Source the synonym check if wanted
  if(tolower(synonym_check_input) == "y") {
    synonym_check()
  } else {
    cat("Skipping synonym Check \n")
  }
  
  merger()
  
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
    duplicate_species = rgbif_similarity(scientific_names)
  } 
  
  ## combine the two columns into one
  if (a == "n") {
    simCheck_filt_sp = jw_similarity_check()
    scientific_names = c(simCheck_filt_sp[ ,"name", drop = TRUE], simCheck_filt_sp[ ,"similarName", drop = TRUE])
    scientific_names = scientific_names[!duplicated(scientific_names)]
  }
  
  
  # Further use rgbif to check the similar names
  #duplicate_species = rgbif_similarity(filtered_species$scientificName)
  
  ## Remove the duplications from the GBIF check
  filtered_species_final = anti_join(filtered_species, duplicate_species, by = "scientificName") #Double check if this is the correct way to do it
  write.csv(filtered_species_final, "outputs/filtered_species_final.csv", row.names = F, fileEncoding = "UTF-8")
  
  # Stop the filter script timer
  end_filter_time = Sys.time()
  ### Calculate the time
  el_time_filter = format_elapsed_time(start_filter_time, end_filter_time)
  # Print message with elapsed time
  message("Script finished in ", el_time_filter)
}

filtering()
