message("Initiating GBIF Species Retrieval")
# Retrieve all species names within the Tracheophyta phylum
## get taxonomykey for tracheophyta
tracheophytaTaxonKey = name_backbone("Tracheophyta")$usageKey

cat("Collecting GBIF species Keys \n")
## Count amount of unique species instances
gbif_species_count = occ_count(
  taxonKey = tracheophytaTaxonKey,
  facet = "speciesKey",
  facetLimit = 20000,
  hasCoordinate = T,
  geometry = combined_WKT
)$speciesKey


cat("Using Keys to collect names. this may take a while \n")
## Use the count to get the species names
gbif_species = sapply(gbif_species_count, FUN = function(x) {
  occ_count(
    speciesKey = x,
    facet = "scientificName",
    hasCoordinate = T,
    geometry = combined_WKT
  )$scientificName
})

## Make the species names into strings
gbif_species = as.character(do.call(c, gbif_species))
gbif_species = as.data.frame(gbif_species)
gbif_species$gbif_species = trimws(gbif_species$gbif_species)

write.csv(gbif_species_count, "outputs/gbif_retrieval_process/gbif_species_count.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(gbif_species, "outputs/gbif_retrieval_process/gbif_species_names.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(gbif_species, "outputs/origin_lists/gbif_species.csv", row.names = F, fileEncoding = "UTF-8")

cat("Retrieved GBIF Species \n")