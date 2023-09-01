merger = function() {
  # choose Scientifically accepted names
  message("------ Merging ABA, AMBIO, GBIF, and GloNAF lists ------")
  
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
  
  aba_ambio = function() {
    present = function() {
      ## combine present lists and remove duplicates
      arctic_present = union(aba_arctic_present, ambio_arctic_present)
      
      ## Run NA and distinct check
      cat("All arctic_present NAs removed:", any(!is.na(arctic_present)), "\n")
      cat("All arctic_present values are distinct: ", any(!duplicated(arctic_present)), "\n")
      
      ## Write out list
      write.csv(arctic_present, "outputs/filtering_process_outputs/combined_aba_ambio_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        arctic_present
      )
    }
    
    absent = function() {
      # Combine AMBIO and ABA lists
      message("------ Merging ABA and AMBIO lists ------")
      
      
      ## combine absent lists and remove duplicates
      arctic_absent = union(aba_arctic_absent, ambio_arctic_absent)
      
      ## Run NA and distinct check
      cat("All arctic_absent NAs removed:", any(!is.na(arctic_absent)), "\n")
      cat("All arctic_absent values are distinct: ", any(!duplicated(arctic_absent)), "\n")
      
      ## Write out list
      write.csv(arctic_absent, "outputs/filtering_process_outputs/combined_aba_ambio_absent.csv",  row.names = F, fileEncoding = "UTF-8")
      
      return(
        arctic_absent
      )  
    }
    
    return(
      arctic_present = present(),
      arctic_absent = absent() 
    )
  }
  
  
  gbif_aba_ambio = function(present, absent) {
    message("------ combine with GBIF list ------")
    
    present = present
    absent = absent
    
    present = function() {
      ## return rows from GBIF that are in the Arctic
      gbif_arctic_present = merge(gbif_species, 
                                  present, 
                                  by = "scientificName"
      )
      
      ## Write out list
      write.csv(gbif_arctic_present, "outputs/filtering_process_outputs/gbif_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      ## combine GBIF present list and arctic present list 
      gbif_aba_ambio_present = union(gbif_arctic_present, present)
      write.csv(gbif_aba_ambio_present, "outputs/filtering_process_outputs/gbif_aba_ambio_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        gbif_aba_ambio_present
      )
    }
    
    absent = function() {
      ## Remove rows from gbif_species data frame that have matching values to the arctic_present species
      gbif_arctic_absent = anti_join(gbif_species, 
                                     present, 
                                     by = "scientificName"
      )
      
      ## Write out list
      write.csv(gbif_arctic_absent, "outputs/filtering_process_outputs/gbif_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      ## Combine GBIF and arctic_absent list and at the same time remove duplicates
      gbif_aba_ambio_absent = union(gbif_arctic_absent, absent)
      
      ## Check how many are similar
      cat("Number of rows in GBIF that are similar to the arctic_absent list: ", nrow(anti_join(gbif_species, present, by = "scientificName")) - nrow(anti_join(gbif_arctic_absent, absent, by = "scientificName")), "\n")
      
      write.csv(gbif_aba_ambio_absent, "outputs/filtering_process_outputs/gbif_aba_ambio_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        gbif_aba_ambio_absent
      )
    }
    
    
    
    cat("gbif_arctic_present and gbif_arctic_absent lists have been made \n")
    
    return(
      gbif_aba_ambio_present = present(),
      gbif_aba_ambio_absent = absent()
    )
  }
  
  
  glonaf_gbif_aba_ambio = function(present, absent) {
    message("------ Merging the combined ABA, AMBIO, and GBIF list with the GloNAF list ------")
    
    present = present
    absent = absent
    
    present = function() {
      ## This will make a new data frame with only the names that exist in both data frames
      ## This means that I will get a list of all the glonaf species that are known to be in the Arctic
      glonaf_arctic_present = merge(glonaf_species, 
                                    present, 
                                    by = "scientificName"
      )
      
      write.csv(glonaf_arctic_present, "outputs/filtering_process_outputs/glonaf_arctic_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      ## Combine glonaf with gbif_aba_ambio arctic present list
      glonaf_gbif_ambio_aba_present = union(glonaf_arctic_present, present)
      write.csv(glonaf_gbif_ambio_aba_present, "outputs/filtering_process_outputs/glonaf_gbif_ambio_aba_present.csv", row.names = F, fileEncoding = "UTF-8")
      
      return(
        glonaf_gbif_ambio_aba_present
      )
    }
    
    absent = function() {
      
      ## Remove arctic present species from the glonaf list to get all species in the glonaf list that are absent in the arctic 
      glonaf_arctic_absent = anti_join(glonaf_species, 
                                       present, 
                                       by = "scientificName"
      )
      
      write.csv(glonaf_arctic_absent, "outputs/filtering_process_outputs/glonaf_arctic_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      cat("All glonaf_arctic_absent NAs removed:", any(!is.na(glonaf_arctic_absent)), "\n")
      cat("All glonaf_arctic_absent are distinct:", any(!duplicated(glonaf_arctic_absent)), "\n")
      
      ## Combine glonaf with gbif_aba_ambio arctic absent list
      glonaf_gbif_ambio_aba_absent = union(glonaf_arctic_absent, absent)
      write.csv(glonaf_gbif_ambio_aba_absent, "outputs/filtering_process_outputs/glonaf_gbif_ambio_aba_absent.csv", row.names = F, fileEncoding = "UTF-8")
      
      
      return(
        glonaf_gbif_ambio_aba_absent
      )
    }
    
    return(
      glonaf_gbif_ambio_aba_present = present(),
      glonaf_gbif_ambio_aba_absent = absent()
    )
    
  }
  
  
  filter = function(absent) {
    message("------ Convert to filtered_species list ------")
    
    filtered_species = glonaf_gbif_ambio_aba_absent
    
    ## Check for NAs and duplications
    cat("All filtered_species NAs removed:", any(!is.na(filtered_species)), "\n")
    cat("All filtered_species are distinct:", any(!duplicated(filtered_species)), "\n")
    
    ## Filter out only Genus names
    # Remove rows where the scientificName column contains only a single word
    filtered_species = filtered_species %>%
      filter(!str_detect(scientificName, "^\\S+$"))
    
    # Make into dataframea again
    filtered_species = data.frame(filtered_species)
    
    # Replace odd multiplication sign with x in order to better display the final list
    filtered_species$scientificName = lapply(filtered_species$scientificName, function(x) gsub("(\\w) × (\\w)", "\\1 x \\2", x))
    
    # Unlist filtered_species to make it a vector
    filtered_species$scientificName = unlist(filtered_species$scientificName)
    
    return(filtered_species)
  }
  
  aba_ambio()
  gbif_aba_ambio(arctic_present, arctic_absent)
  glonaf_gbif_aba_ambio(gbif_aba_ambio_present, gbif_aba_ambio_absent)
  
  
  return(
    filtered_species = filter(glonaf_gbif_ambio_aba_absent)
  )
}
