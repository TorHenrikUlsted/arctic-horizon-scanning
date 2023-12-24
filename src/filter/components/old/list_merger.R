merger = function() {
    
    gbif_aba_ambio_p = function(arctic_present) {
      message("------ Combine with GBIF list ------")
      
      gbif_species = results$gbif_species
      
      cat(cc$paleTurquoise("Creating GBIF Arctic present species \n"))
      gbif_p = merge(gbif_species, arctic_present, by = "scientificName")
      
      cat(yellow("writing out csv file to: \n", "./outputs/filtering/filtering_process_outputs/gbif_arctic_present.csv \n"))
      write.csv(gbif_p, "./outputs/filtering/filtering_process_outputs/gbif_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      cat("Combining the Arctic Present list with the GBIF present list \n")
      gbif_aba_ambio_p = dplyr::union(gbif_p, arctic_present)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(gbif_aba_ambio_p)) == T) {green(any(!is.na(gbif_aba_ambio_p)))} else {red(any(!is.na(gbif_aba_ambio_p)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(gbif_aba_ambio_p)) == T) {green(any(!duplicated(gbif_aba_ambio_p)))} else {red(any(!duplicated(gbif_aba_ambio_p)))}, "\n")
      cat("Duplicated names found: ", green(nrow(gbif_p) + nrow(arctic_present) - nrow(gbif_aba_ambio_p)), "\n")
      
      cat(yellow("Writing out list to: \n", "./outputs/filtering/filtering_process_outputs/gbif_aba_ambio_present.csv \n\n"))
      write.csv(gbif_aba_ambio_p, "./outputs/filtering/filtering_process_outputs/gbif_aba_ambio_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        gbif_aba_ambio_p
      )
    }
    
    gbif_aba_ambio_a = function(arctic_present, arctic_absent) {
      gbif_species = results$gbif_species
      
      cat(cc$paleTurquoise("Creating GBIF Arctic absent species list \n"))
      ## Remove rows from GBIF data frame that have matching values to the arctic_present species
      gbif_a = anti_join(gbif_species, arctic_present, by = "scientificName")
      
      cat(yellow("Writing out list to: \n", "./outputs/filtering/filtering_process_outputs/gbif_arctic_absent.csv \n"))
      write.csv(gbif_a, "./outputs/filtering/filtering_process_outputs/gbif_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")

      cat("Combining GBIF absent list with Arctic absent list ")
      gbif_aba_ambio_a = dplyr::union(gbif_a, arctic_absent)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(gbif_aba_ambio_a)) == T) {green(any(!is.na(gbif_aba_ambio_a)))} else {red(any(!is.na(gbif_aba_ambio_a)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(gbif_aba_ambio_a)) == T) {green(any(!duplicated(gbif_aba_ambio_a)))} else {red(any(!duplicated(gbif_aba_ambio_a)))}, "\n")
      cat("Duplicated names found: ", green(nrow(gbif_a) + nrow(arctic_absent) - nrow(gbif_aba_ambio_a)), "\n")
      
      cat(yellow("Writing out csv file to: \n", "./outputs/filtering/filtering_process_outputs/gbif_aba_ambio_absent.csv \n\n"))
      write.csv(gbif_aba_ambio_a, "./outputs/filtering/filtering_process_outputs/gbif_aba_ambio_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        gbif_aba_ambio_a
      )
    }
  
    
    glonaf_gbif_aba_ambio_p = function(arctic_present) {
      message("------ Merging the combined ABA, AMBIO, and GBIF list with the GloNAF list ------")
      glonaf_species = results$glonaf_species
      
      cat(cc$paleTurquoise("Creating GloNAF Arctic present species list \n"))
      ## This will make a new data frame with only the names that exist in both data frames
      ## This means that I will get a list of all the glonaf species that are known to be in the Arctic
      glonaf_p = merge(glonaf_species, arctic_present, by = "scientificName")
      
      cat(yellow("Writing out csv file to: \n", "./outputs/filtering/filtering_process_outputs/glonaf_arctic_present.csv \n"))
      write.csv(glonaf_p, "./outputs/filtering/filtering_process_outputs/glonaf_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      ## Combine glonaf with gbif_aba_ambio arctic present list
      cat("Creating a combined GloNAF, GBIF, ABA, AMBIO arctic present list \n")
      glonaf_gbif_ambio_aba_p = dplyr::union(glonaf_p, arctic_present)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(glonaf_gbif_ambio_aba_p)) == T) {green(any(!is.na(glonaf_gbif_ambio_aba_p)))} else {red(any(!is.na(glonaf_gbif_ambio_aba_p)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(glonaf_gbif_ambio_aba_p)) == T) {green(any(!duplicated(glonaf_gbif_ambio_aba_p)))} else {red(any(!duplicated(glonaf_gbif_ambio_aba_p)))}, "\n")
      cat("Duplicated names found: ", green(nrow(glonaf_p) + nrow(arctic_present) - nrow(glonaf_gbif_ambio_aba_p)), "\n")
      
      cat(yellow("Writing out csv file to: \n", "./outputs/filtering/filtering_process_outputs/glonaf_gbif_ambio_aba_present.csv \n\n"))
      write.csv(glonaf_gbif_ambio_aba_p, "./outputs/filtering/filtering_process_outputs/glonaf_gbif_ambio_aba_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        glonaf_gbif_ambio_aba_p
      )
    }
    
    glonaf_gbif_aba_ambio_a = function(arctic_present, arctic_absent) {
      glonaf_species = results$glonaf_species
      
      cat(cc$paleTurquoise("Creating GloNAF Arctic absent species list \n"))
      ## Remove arctic present species from the glonaf list to get all species in the glonaf list that are absent in the arctic 
      glonaf_a = anti_join(glonaf_species, arctic_present, by = "scientificName")
      
      cat(yellow("Writing csv file to: \n", "./outputs/filtering/filtering_process_outputs/glonaf_arctic_absent.csv \n"))
      write.csv(glonaf_a, "./outputs/filtering/filtering_process_outputs/glonaf_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      ## Combine glonaf with gbif_aba_ambio arctic absent list
      cat("Combining GloNAF absent list with the GBIF, ABA, AMBIO combined list \n")
      glonaf_gbif_ambio_aba_a = dplyr::union(glonaf_a, arctic_absent)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(glonaf_gbif_ambio_aba_a)) == T) {green(any(!is.na(glonaf_gbif_ambio_aba_a)))} else {red(any(!is.na(glonaf_gbif_ambio_aba_a)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(glonaf_gbif_ambio_aba_a)) == T) {green(any(!duplicated(glonaf_gbif_ambio_aba_a)))} else {red(any(!duplicated(glonaf_gbif_ambio_aba_a)))}, "\n")
      cat("Duplicated names found: ", green(nrow(glonaf_a) + nrow(arctic_absent) - nrow(glonaf_gbif_ambio_aba_a)), "\n")
      
      cat(yellow("Writing out csv file to: \n", "./outputs/filtering/filtering_process_outputs/glonaf_gbif_ambio_aba_absent.csv \n\n"))
      write.csv(glonaf_gbif_ambio_aba_a, "./outputs/filtering/filtering_process_outputs/glonaf_gbif_ambio_aba_absent.csv", row.names = F, fileEncoding = "UTF-8")
      return(
        glonaf_gbif_ambio_aba_a
      )
    }
    
  
  
  filter_species = function(arctic_absent) {
    message("------ Convert to filtered_species list ------")
    
    filtered_species = arctic_absent
    
    ## Check for NAs and duplications
    ## Run NA and distinct check
    cat("All filtered_species NAs removed:", if(any(!is.na(filtered_species)) == T) {green(any(!is.na(filtered_species)))} else {red(any(!is.na(filtered_species)))}, "\n")
    cat("All filtered_species values are distinct: ", if(any(!duplicated(filtered_species)) == T) {green(any(!duplicated(filtered_species)))} else {red(any(!duplicated(filtered_species)))}, "\n")
    
    ## Filter out only Genus names
    # Remove rows where the scientificName column contains only a single word
    cat("Filtering out only Genus names \n")
    filtered_species = filtered_species %>%
      filter(!str_detect(scientificName, "^\\S+$"))
    
    # Make into dataframea again
    cat("Creating data frame and changing multiplication sign to letter x \n")
    filtered_species = data.frame(filtered_species)
    
    # Replace odd multiplication sign with x in order to better display the final list
    filtered_species$scientificName = lapply(filtered_species$scientificName, function(x) gsub("(\\w) × (\\w)", "\\1 x \\2", x))
    
    ## Run NA and distinct check
    cat("All filtered_species NAs removed:", if(any(!is.na(filtered_species)) == T) {green(any(!is.na(filtered_species)))} else {red(any(!is.na(filtered_species)))}, "\n")
    cat("All filtered_species values are distinct: ", if(any(!duplicated(filtered_species)) == T) {green(any(!duplicated(filtered_species)))} else {red(any(!duplicated(filtered_species)))}, "\n")
    
    # Unlist filtered_species to make it a vector
    cat("Unlisting data frame for further use \n")
    filtered_species$scientificName = unlist(filtered_species$scientificName)
    filtered_species = na.omit(filtered_species)
    
    cat(yellow("Writing out csv file to: \n", "./outputs/filtering/filtering_process_outputs/filtered_species.csv \n\n"))
    write.csv(filtered_species, "./outputs/filtering/filtering_process_outputs/filtered_species.csv", row.names = F, fileEncoding = "UTF-8")
    
    cat(cc$aquamarine("Scroll up and check if all boolean operators are green / TRUE. If not, something is not working correctly and debugging is necessary. \n\n"))
    
    return(filtered_species)
  }
  
  aba_ambio_p = aba_ambio_p()
  aba_ambio_a = aba_ambio_a()
  
  gbif_aba_ambio_p = gbif_aba_ambio_p(aba_ambio_p)
  gbif_aba_ambio_a = gbif_aba_ambio_a(aba_ambio_p, aba_ambio_a)
  
  glonaf_gbif_aba_ambio_p = glonaf_gbif_aba_ambio_p(gbif_aba_ambio_p)
  glonaf_gbif_aba_ambio_a = glonaf_gbif_aba_ambio_a(gbif_aba_ambio_p, gbif_aba_ambio_a)
  
  return(
    filtered_species = filter_species(glonaf_gbif_aba_ambio_a)
  )
}
