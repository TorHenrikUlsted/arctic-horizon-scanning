merger = function() {
  results = list()
  cc = custom_colors()
  
  # choose Scientifically accepted names
  message("------ Merging all lists to create the filtered_species list ------")
  
  ## Define components and their associated file paths
  components = list(
    aba = list(
      present = "outputs/wfo_one_outputs/wfo_one_aba_arctic_present.csv",
      absent = "outputs/wfo_one_outputs/wfo_one_aba_arctic_absent.csv"
    ),
    ambio = list(
      present = "outputs/wfo_one_outputs/wfo_one_ambio_arctic_present.csv",
      absent = "outputs/wfo_one_outputs/wfo_one_ambio_arctic_absent.csv"
    ),
    gbif = list(species = "outputs/wfo_one_outputs/wfo_one_gbif_species.csv"),
    glonaf = list(species = "outputs/wfo_one_outputs/wfo_one_glonaf_species.csv")
  )
  
  ## Iterate over components
  for (component in names(components)) {
    if (component != "gbif" && component != "glonaf") {
      present_data = read.csv(components[[component]]$present)
      absent_data = read.csv(components[[component]]$absent)
      
      ## Select the scientific names and remove NA
      present_data = present_data %>% 
        select(scientificName) %>% 
        filter(!is.na(scientificName))
      cat(paste0("All ", component, "_arctic_present NAs removed: ", if(any(!is.na(present_data)) == T) {green(any(!is.na(present_data)))} else {red(any(!is.na(present_data)))}, "\n"))
      
      ## Remove duplicates
      present_data = distinct(present_data)
      cat(paste0("All ", component, "_arctic_present are distinct: ", if(any(!duplicated(present_data)) == T) {green(any(!duplicated(present_data)))} else {red(any(!duplicated(present_data)))}, "\n"))
      
      ## Select the scientific names and remove NA
      absent_data = absent_data %>% 
        select(scientificName) %>% 
        filter(!is.na(scientificName))
      cat(paste0("All ", component, "_arctic_absent NAs removed: ", if(any(!is.na(absent_data)) == T) {green(any(!is.na(absent_data)))} else {red(any(!is.na(absent_data)))}, "\n"))
      
      ## Remove duplicates
      absent_data = distinct(absent_data)
      cat(paste0("All ", component, "_arctic_absent are distinct: ", if(any(!duplicated(absent_data)) == T) {green(any(!duplicated(absent_data)))} else {red(any(!duplicated(absent_data)))}, "\n"))
      
      # store the results in a list
      results[[paste0(component, "_arctic_present")]] = present_data
      results[[paste0(component, "_arctic_absent")]] = absent_data
    } else {
      species_data = read.csv(components[[component]]$species)
      species_data = species_data %>% 
        select(scientificName) %>% 
        filter(!is.na(scientificName))
      
      ## Remove duplicates
      species_data = distinct(species_data)
      
      results[[paste0(component, "_species")]] = species_data
    }
  }
  cat(cc$lightSteelBlue("\n The neccessary lists are: aba_arctic_present/absent, ambio_arctic_present/absent, gbif_species, glonaf_spcies \n"))
  cat(cc$lightSteelBlue("These are the registered species lists: "), "\n-", paste(names(results), sep = " ", collapse="\n- "), "\n\n")
  
  # Combine ABA and AMBIO data frames
  ## No need for parameters as R allows for lexical scoping -- meaning access to parental variables
    aba_ambio_p = function() {
      message("------ Merging ABA and AMBIO lists ------")
      
      aba_arctic_present = results$aba_arctic_present
      ambio_arctic_present  = results$ambio_arctic_present
      
      cat(cc$paleTurquoise("Combining ABA and AMBIO present lists \n"))
      ## combine present lists and remove duplicates
      arctic_present = dplyr::union(aba_arctic_present, ambio_arctic_present)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(arctic_present)) == T) {green(any(!is.na(arctic_present)))} else {red(any(!is.na(arctic_present)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(arctic_present)) == T) {green(any(!duplicated(arctic_present)))} else {red(any(!duplicated(arctic_present)))}, "\n")
      cat("Duplicated names found: ", green(nrow(aba_arctic_present) + nrow(ambio_arctic_present) - nrow(arctic_present)), "\n")

      cat(yellow("writing out csv file to: \n", "outputs/filtering_process_outputs/combined_aba_ambio_present.csv \n\n"))
      write.csv(arctic_present, "outputs/filtering_process_outputs/combined_aba_ambio_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        arctic_present
      )
    }
    
    aba_ambio_a = function() {
      aba_arctic_absent = results$aba_arctic_absent
      ambio_arctic_absent  = results$ambio_arctic_absent
      
      cat(cc$paleTurquoise("Combining ABA and AMBIO absent lists \n"))
      ## combine absent lists and remove duplicates
      arctic_absent = dplyr::union(aba_arctic_absent, ambio_arctic_absent)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(arctic_absent)) == T) {green(any(!is.na(arctic_absent)))} else {red(any(!is.na(arctic_absent)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(arctic_absent)) == T) {green(any(!duplicated(arctic_absent)))} else {red(any(!duplicated(arctic_absent)))}, "\n")
      cat("Duplicated names found: ", green(nrow(aba_arctic_absent) + nrow(ambio_arctic_absent) - nrow(arctic_absent)), "\n")
      
      ## Write out list
      cat(yellow("writing out csv file to: \n", "outputs/filtering_process_outputs/combined_aba_ambio_absent.csv \n\n"))
      write.csv(arctic_absent, "outputs/filtering_process_outputs/combined_aba_ambio_absent.csv",  row.names = F, fileEncoding = "UTF-8")
      
      return(
        arctic_absent
      )  
    }
   
    
    gbif_aba_ambio_p = function(arctic_present) {
      message("------ Combine with GBIF list ------")
      
      gbif_species = results$gbif_species
      
      cat(cc$paleTurquoise("Creating GBIF Arctic present species \n"))
      gbif_p = merge(gbif_species, arctic_present, by = "scientificName")
      
      cat(yellow("writing out csv file to: \n", "outputs/filtering_process_outputs/gbif_arctic_present.csv \n"))
      write.csv(gbif_p, "outputs/filtering_process_outputs/gbif_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      cat("Combining the Arctic Present list with the GBIF present list \n")
      gbif_aba_ambio_p = dplyr::union(gbif_p, arctic_present)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(gbif_aba_ambio_p)) == T) {green(any(!is.na(gbif_aba_ambio_p)))} else {red(any(!is.na(gbif_aba_ambio_p)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(gbif_aba_ambio_p)) == T) {green(any(!duplicated(gbif_aba_ambio_p)))} else {red(any(!duplicated(gbif_aba_ambio_p)))}, "\n")
      cat("Duplicated names found: ", green(nrow(gbif_p) + nrow(arctic_present) - nrow(gbif_aba_ambio_p)), "\n")
      
      cat(yellow("Writing out list to: \n", "outputs/filtering_process_outputs/gbif_aba_ambio_present.csv \n\n"))
      write.csv(gbif_aba_ambio_p, "outputs/filtering_process_outputs/gbif_aba_ambio_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        gbif_aba_ambio_p
      )
    }
    
    gbif_aba_ambio_a = function(arctic_present, arctic_absent) {
      gbif_species = results$gbif_species
      
      cat(cc$paleTurquoise("Creating GBIF Arctic absent species list \n"))
      ## Remove rows from GBIF data frame that have matching values to the arctic_present species
      gbif_a = anti_join(gbif_species, arctic_present, by = "scientificName")
      
      cat(yellow("Writing out list to: \n", "outputs/filtering_process_outputs/gbif_arctic_absent.csv \n"))
      write.csv(gbif_a, "outputs/filtering_process_outputs/gbif_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")

      cat("Combining GBIF absent list with Arctic absent list ")
      gbif_aba_ambio_a = dplyr::union(gbif_a, arctic_absent)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(gbif_aba_ambio_a)) == T) {green(any(!is.na(gbif_aba_ambio_a)))} else {red(any(!is.na(gbif_aba_ambio_a)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(gbif_aba_ambio_a)) == T) {green(any(!duplicated(gbif_aba_ambio_a)))} else {red(any(!duplicated(gbif_aba_ambio_a)))}, "\n")
      cat("Duplicated names found: ", green(nrow(gbif_a) + nrow(arctic_absent) - nrow(gbif_aba_ambio_a)), "\n")
      
      cat(yellow("Writing out csv file to: \n", "outputs/filtering_process_outputs/gbif_aba_ambio_absent.csv \n\n"))
      write.csv(gbif_aba_ambio_a, "outputs/filtering_process_outputs/gbif_aba_ambio_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
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
      
      cat(yellow("Writing out csv file to: \n", "outputs/filtering_process_outputs/glonaf_arctic_present.csv \n"))
      write.csv(glonaf_p, "outputs/filtering_process_outputs/glonaf_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      ## Combine glonaf with gbif_aba_ambio arctic present list
      cat("Creating a combined GloNAF, GBIF, ABA, AMBIO arctic present list \n")
      glonaf_gbif_ambio_aba_p = dplyr::union(glonaf_p, arctic_present)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(glonaf_gbif_ambio_aba_p)) == T) {green(any(!is.na(glonaf_gbif_ambio_aba_p)))} else {red(any(!is.na(glonaf_gbif_ambio_aba_p)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(glonaf_gbif_ambio_aba_p)) == T) {green(any(!duplicated(glonaf_gbif_ambio_aba_p)))} else {red(any(!duplicated(glonaf_gbif_ambio_aba_p)))}, "\n")
      cat("Duplicated names found: ", green(nrow(glonaf_p) + nrow(arctic_present) - nrow(glonaf_gbif_ambio_aba_p)), "\n")
      
      cat(yellow("Writing out csv file to: \n", "outputs/filtering_process_outputs/glonaf_gbif_ambio_aba_present.csv \n\n"))
      write.csv(glonaf_gbif_ambio_aba_p, "outputs/filtering_process_outputs/glonaf_gbif_ambio_aba_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        glonaf_gbif_ambio_aba_p
      )
    }
    
    glonaf_gbif_aba_ambio_a = function(arctic_present, arctic_absent) {
      glonaf_species = results$glonaf_species
      
      cat(cc$paleTurquoise("Creating GloNAF Arctic absent species list \n"))
      ## Remove arctic present species from the glonaf list to get all species in the glonaf list that are absent in the arctic 
      glonaf_a = anti_join(glonaf_species, arctic_present, by = "scientificName")
      
      cat(yellow("Writing csv file to: \n", "outputs/filtering_process_outputs/glonaf_arctic_absent.csv \n"))
      write.csv(glonaf_a, "outputs/filtering_process_outputs/glonaf_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      ## Combine glonaf with gbif_aba_ambio arctic absent list
      cat("Combining GloNAF absent list with the GBIF, ABA, AMBIO combined list \n")
      glonaf_gbif_ambio_aba_a = dplyr::union(glonaf_a, arctic_absent)
      
      ## Run NA and distinct check
      cat("All NAs removed:", if(any(!is.na(glonaf_gbif_ambio_aba_a)) == T) {green(any(!is.na(glonaf_gbif_ambio_aba_a)))} else {red(any(!is.na(glonaf_gbif_ambio_aba_a)))}, "\n")
      cat("All values are distinct: ", if(any(!duplicated(glonaf_gbif_ambio_aba_a)) == T) {green(any(!duplicated(glonaf_gbif_ambio_aba_a)))} else {red(any(!duplicated(glonaf_gbif_ambio_aba_a)))}, "\n")
      cat("Duplicated names found: ", green(nrow(glonaf_a) + nrow(arctic_absent) - nrow(glonaf_gbif_ambio_aba_a)), "\n")
      
      cat(yellow("Writing out csv file to: \n", "outputs/filtering_process_outputs/glonaf_gbif_ambio_aba_absent.csv \n\n"))
      write.csv(glonaf_gbif_ambio_aba_a, "outputs/filtering_process_outputs/glonaf_gbif_ambio_aba_absent.csv", row.names = F, fileEncoding = "UTF-8")
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
    
    cat(yellow("Writing out csv file to: \n", "outputs/filtering_process_outputs/filtered_species.csv \n\n"))
    write.csv(filtered_species, "outputs/filtering_process_outputs/filtered_species.csv", row.names = F, fileEncoding = "UTF-8")
    
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
