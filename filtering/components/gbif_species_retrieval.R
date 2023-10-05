sp_retrieve = function() {  
  message("Initiating GBIF Species Retrieval")
  # Retrieve all species names within the Tracheophyta phylum
  core_usage = detectCores()/4
  cc = custom_colors()
  
  ## get taxonomykey for tracheophyta
  tracheophytaTaxonKey = name_backbone("Tracheophyta")$usageKey
  
  cat("Collecting GBIF species Keys \n")
  ## Count amount of unique species instances
  gbif_species_count = occ_count(
    taxonKey = tracheophytaTaxonKey,
    facet = "speciesKey",
    facetLimit = 100000,
    hasCoordinate = T,
    geometry = combined_WKT
  )$speciesKey
  
  cat(yellow("Writing out gbif_species_count to: \n", "./outputs/filtering/gbif_retrieval_process/gbif_species_count.csv \n"))
  cat("All species names unique: ", if (any(!duplicated(gbif_species_count)) == T) {green(any(!duplicated(gbif_species_count)))} else {red(any(!duplicated(gbif_species_count)))}, "\n")
  write.csv(gbif_species_count, "./outputs/filtering/gbif_retrieval_process/gbif_species_count.csv", row.names = F, fileEncoding = "UTF-8")
  
  
  cat("Using Keys to collect names. this may take a while \n")
  ## Use the count to get the species names
  # Set up a progress handler
  handlers("txtprogressbar")
  
  # Use with_progress instead of directly calling mclapply
  with_progress({
    # Set up a progress bar
    p = progressor(steps = length(gbif_species_count))
    
  gbif_species = mclapply(gbif_species_count, FUN = function(x) {
  # Update progress bar
    p()
    
    occ_count(
      speciesKey = x,
      facet = "scientificName",
      hasCoordinate = T,
      geometry = combined_WKT
    )$scientificName
  }, mc.cores = core_usage)
  })
  
  ## Make the species names into strings
  gbif_species = as.character(do.call(c, gbif_species))
  gbif_species = data.frame(gbif_species = trimws(gbif_species))
  
  cat(yellow("Writing out gbif_species to: \n", "./outputs/filtering/gbif_retrieval_process/gbif_species_names.csv ...and \n", "./outputs/filtering/origin_lists/gbif_species.csv \n"))
  write.csv(gbif_species, "./outputs/filtering/gbif_retrieval_process/gbif_species_names.csv", row.names = F, fileEncoding = "UTF-8")
  write.csv(gbif_species, "./outputs/filtering/origin_lists/gbif_species.csv", row.names = F, fileEncoding = "UTF-8")
  
  cat("Starting GBIF similarity check \n")
  scientific_names = unlist(gbif_species$gbif_species)
  gbif_species_check = rgbif_simliarity_check(scientific_names)
  cat("Writing out similarity check csv to: \n", "./outputs/filtering/similarity_check_outputs/gbif_species_check.csv \n")
  write.csv(gbif_species_check, "./outputs/filtering/similarity_check_outputs/gbif_species_check.csv", row.names = F, fileEncoding = "UTF-8")
  
  a = ""
  cat(cc$aquamarine("If you recently updated shape files or similar that will affect GBIF geometry, then you need to download a new occurence file from GBIF. \n", "or else just download the existing gbif list. \n"))
  
  while(a != "y" && a != "n") {
    a = readline("Do you want to download a new occurence file [y], or import from csv? \n")
    
    if (a != "y" && a != "n") {
      cat("Invalid response. Please enter 'y' or 'n'.\n")
    }
  }
  ## Use the filtered species list with all the names
  if (a == "y") {
    source("./filtering/components/gbif_occ_retrieval.R")
    sp_occ_df = occ_retrieval(gbif_species_check);
    
    source("./filtering/components/gbif_sp_crop.R")
    sp_occ_cropped = crop_species(sp_occ_df);
  } 
  
  if (a == "n") {
    cat("Loading occurence data from file \n")
    sp_cavm = read.csv("./outputs/filtering/gbif_retrieval_process/gbif_boreal_species.csv")
    sp_boreal = read.csv("./outputs/filtering/gbif_retrieval_process/gbif_boreal_species.csv")
    
    sp_occ_cropped = list(sp_cavm, sp_boreal)
  }
  
  return(sp_occ_cropped)
}