# here will the filtering take place by source("components/filename.r")
# here i will filter all GBIF species from arctic and boreal against ABA and AMBIO
# Also after that i will filter them against the glonaf list in the end, and all glonaf species that are not in the GBIF mix will be removed from the project.




# conduct a synonym check for all datasets
library(WorldFlora)
# Conduct a synonym check using WFO
glonaf_species_wfo = WFO.match(spec.data = glonaf_species, spec.name = "standardized_name", WFO.file = "resources/wfo_classification.csv", verbose = T, counter = 100)
message("WFO completed the match")