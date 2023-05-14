library(rgbif)
source("components/regions_data_import.R")

# Retrieve all species names within the Tracheophyta phylum
## Download GBIF Backbone Taxonomy
cat("Download GBIF Backbone Taxonomy here: ", "https://doi.org/10.15468/39omei", "\n")
## Read backbone
gbif_backbone = read.delim("resources/Taxon.tsv", stringsAsFactors = F)
cat("Completed reading GBIF Backbone Checklist")

## Extract rows for species rank
gbif_species = subset(gbif_backbone, taxonRank == "species")

## Extract rows for Tracheophyta phylum
gbif_tracheophyta_species = subset(gbif_species, phylum == "Tracheophyta")

## Extract the canonical names
gbif_tracheophyta_species = select(gbif_tracheophyta_species, Species = scientificName)


# Run WFO synonym check


## ----------------- Comment out start ------------------ ##

# # Download the full dataset
# taxonKey = name_backbone(name = "Tracheophyta")$usageKey
# ## Test data download
# occ_download_prep(pred("taxonKey", taxonKey), 
#                   pred("hasGeospatialIssue", FALSE),
#                   pred("hasCoordinate", TRUE),
#                   pred("geometry", combined_WKT),
#                   pred("occurrenceStatus","PRESENT"),
#                   pred_in("basisOfRecord", c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "OCCURENCE")),
#                   format = "SIMPLE_CSV"
# )

## Download the data 
# gbif_species = occ_download(pred("taxonKey", taxonKey), 
#                                pred("hasGeospatialIssue", FALSE),
#                                pred("hasCoordinate", TRUE),
#                                pred("geometry", combined_WKT),
#                                pred("occurrenceStatus","PRESENT"),
#                                pred_in("basisOfRecord", c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "OCCURENCE")),
#                                format = "SIMPLE_CSV"
# )

## check status
# occ_download_wait(gbif_species)
# 
# ## get the download Data and import to create dataframe
# gbif_species_file = occ_download_get(gbif_species, path = "resources", overwrite = T)
# gbif_species_df = occ_download_import(gbif_species, path = "resources" )

## USE THIS IF ALREADY DOWNLOADED
#gbif_species_df = occ_download_import(as.download("resources/0083144-230224095556074.zip"))

## ----------------- Comment out end ------------------ ##