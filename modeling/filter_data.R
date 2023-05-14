# here will the filtering take place by source("components/filename.r")
# here i will filter all GBIF species from arctic and boreal against ABA and AMBIO
# Also after that i will filter them against the glonaf list in the end, and all GBIF species not in the glonaf list will be removed?

## Source the wrangled files, GBIF will have to be run separately because of download.
source("components/aba_wrangling.R")
source("components/ambio_wrangling.R")
message("sourcing ABA and AMBIO completed")

# Combine ambio and aba lists
## This process will not remove all duplicate names
## merge absent lists
arctic_absent = aba_arctic_absent %>% 
  full_join(ambio_arctic_absent, by = "Species_SubSpecies") %>% 
  distinct(Species_SubSpecies, .keep_all = TRUE)

## merge present lists
arctic_present = aba_arctic_present %>% 
  full_join(ambio_arctic_present, by = "Species_SubSpecies") %>% 
  distinct(Species_SubSpecies, .keep_all = TRUE)


## Remember to have the GBIF occ_download() function commented out before running it here
source("components/gbif_data_retreival.R")


# conduct a synonym check for all datasets
library(WorldFlora)
# Conduct a synonym check using WFO
glonaf_species_wfo = WFO.match(spec.data = glonaf_species, spec.name = "standardized_name", WFO.file = "resources/wfo_classification.csv", verbose = T, counter = 100)
message("WFO completed the match")
?WFO.match()
