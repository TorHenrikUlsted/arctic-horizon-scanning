message("Initiating GBIF Occurrence retrieval")
 # Download the full dataset
 taxonKey = name_backbone(name = "Tracheophyta")$usageKey
 ## Test data download
 occ_download_prep(pred("taxonKey", taxonKey), 
                   pred("hasGeospatialIssue", FALSE),
                   pred("hasCoordinate", TRUE),
                   pred("geometry", combined_WKT),
                   pred("occurrenceStatus","PRESENT"),
                   pred_in("basisOfRecord", c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "OCCURENCE")),
                   format = "SIMPLE_CSV"
 )

## Download the data 
 gbif_species = occ_download(pred("taxonKey", taxonKey), 
                                pred("hasGeospatialIssue", FALSE),
                                pred("hasCoordinate", TRUE),
                                pred("geometry", combined_WKT),
                                pred("occurrenceStatus","PRESENT"),
                                pred_in("basisOfRecord", c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "OCCURENCE")),
                                format = "SIMPLE_CSV"
 )

## check status
occ_download_wait(gbif_species)

## get the download Data and import to create dataframe
gbif_species_file = occ_download_get(gbif_species, path = "resources", overwrite = T)
gbif_species_df = occ_download_import(gbif_species, path = "resources" )

## USE THIS IF ALREADY DOWNLOADED
gbif_species_df = occ_download_import(as.download("resources/0083144-230224095556074.zip"))
cat("Retrieved GBIF occurrences \n")