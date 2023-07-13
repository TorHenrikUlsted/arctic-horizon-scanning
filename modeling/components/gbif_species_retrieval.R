message("Initiating GBIF Species Retrieval")
# Retrieve all species names within the Tracheophyta phylum
## get taxonomykey for tracheophyta
tracheophytaTaxonKey = name_backbone("Tracheophyta")$usageKey

## Count amount of unique species instances
gbif_species_count = occ_count(
  taxonKey = tracheophytaTaxonKey,
  facet = "speciesKey",
  facetLimit = 10000,
  hasCoordinate = T,
  geometry = combined_WKT
)$speciesKey

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

cat("Retrieved GBIF Species \n")