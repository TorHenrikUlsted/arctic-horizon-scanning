get_gbif_keys_region <- function(region, region.name, taxon = "Tracheophyta", verbose = F) {
  cat(blue("Initiating GBIF key by region acquisition protocol. \n"))
  ## get taxonomykey for tracheophyta
  taxon_key <- name_backbone(taxon)$usageKey
  
  cat("Collecting GBIF species Keys \n")
  ## Count amount of unique species instances
  gbif_species_count <- occ_count(
    taxonKey = taxon_key,
    facet = "speciesKey",
    facetLimit = 100000,
    hasCoordinate = T,
    geometry = region
  )$speciesKey
  
  if (any(duplicated(gbif_species_count)) == T) {
    cat(cc$lightCoral("Some keys are duplicated. \n"))
    cat("Attempting to remove... \n")
    
    gbif_species_count <- unique(gbif_species_count)
    
    if (any(duplicated(gbif_species_count)) == T) cat(cc$lightCoral("failed to remove duplications. \n")) else cat(cc$lightGreen("successfully removed duplications. \n"))
    
  } else {
    cat(cc$lightGreen("All keys are unique. \n"))
  }
  
  cat("Number of species keys found", cc$lightSteelBlue(length(gbif_species_count)), "\n")
  
  gbif_count_df <- data.frame(speciesKeys = gbif_species_count)
  print(head(gbif_species_count, 10))
  
  count_savefile <- paste0("./outputs/filter/gbif-acquisition/sp-count-", region.name, ".csv")
  
  cat(yellow("Writing out gbif_species_count to: \n", count_savefile, "\n"))
  create_dir_if("./outputs/filter/gbif-acquisition")
  gbif_count_df <- set_df_utf8(gbif_count_df)
  fwrite(gbif_count_df, count_savefile, row.names = F, bom = T)
  
  cat(cc$lightGreen("GBIF key by region acquisition protocol completed successfully. \n"))
  return(gbif_species_count)
}