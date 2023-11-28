get_gbif_counts <- function(region = NULL, df = NULL) {
  cat("Collecting GBIF species Keys \n")
  
  names(df)[1] <- "scientificName" 
  
  cat("Getting count of species using the wkt: \n", region, "\n")
  
  gbif_counts <- lapply(df$scientificName, function(taxon_name) {
    taxon_key <- name_backbone(taxon_name)$usageKey
    
    gbif_count <- occ_count(
      taxonKey = taxon_key,
      facet = "speciesKey",
      facetLimit = 100000,
      hasCoordinate = T,
      geometry = region
    )$speciesKey
    
    cat(cc$aquamarine(length(unique(gbif_count))), green("unique species keys collected for ", taxon_name, "\n"))
    
    return(gbif_count)
  })
  
  names(gbif_counts) <- df$scientificNames
  
  return(gbif_counts)
}
