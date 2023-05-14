# here will the filtering take place by source("components/filename.r")
# here i will filter all GBIF species from arctic and boreal against ABA and AMBIO
# Also after that i will filter them against the glonaf list in the end, and all GBIF species not in the glonaf list will be removed?

# Source scripts
## Source utility script
source("components/utils.R")
## Source the wrangled ABA and AMBIO lists
abaAmbioList_response = readline("Do you have the abaAmbio_arctic_present/absent CSV files already (y/n)? \n")
source("components/aba_wrangling.R")
source("components/ambio_wrangling.R")
cat("sourcing ABA and AMBIO completed \n")

# Combine ambio and aba lists

## merge present lists
arctic_present = aba_arctic_present %>% 
  full_join(ambio_arctic_present, by = "Species_SubSpecies") %>% 
  distinct(Species_SubSpecies, .keep_all = TRUE)

## merge absent lists
arctic_absent = aba_arctic_absent %>% 
  full_join(ambio_arctic_absent, by = "Species_SubSpecies") %>% 
  distinct(Species_SubSpecies, .keep_all = TRUE)

## write into CSV file
if (abaAmbioList_response == "y") {
  cat("Skipping the creation of abaAmbio CSV files \n")
} else {
  cat("Creating abaAmbio CSV files \n")
  write.csv(arctic_present, "outputs/abaAmbio_arctic_present.csv")
  write.csv(arctic_absent, "outputs/abaAmbio_arctic_absent.csv")
}

## Source the GBIF
## Remember to have the GBIF occ_download() function commented out before running it here
source("components/gbif_data_retrieval.R")
cat("sourcing GBIF data retrieval complete \n")

## Remove rows from gbif_tracheophyta_species data frame that have matching values to the arctic_present species
arctic_present_woAlienInfo = gbif_tracheophyta_species %>% 
  anti_join(arctic_present, by = c("Species" = "Species_SubSpecies"))

## Remove rows from gbif_tracheophyta_species data frame that have matching values to the arctic_absent species
arctic_absent_woAlienInfo = gbif_tracheophyta_species %>% 
  anti_join(arctic_absent, by = c("Species" = "Species_SubSpecies"))





