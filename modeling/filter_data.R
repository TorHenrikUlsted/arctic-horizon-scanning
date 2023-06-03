# Filter all species lists
### Start filter script timer
start_filterScript_time = Sys.time()
## Source scripts ETC 10 min without the synonym check
startTime = Sys.time()
cat("Starting to source scripts. \n")
## Source utility script
source("components/utils.R")
## create the possible commands
simpleListNames = c("aba", "ambio", "gbif", "glonaf")
multiListNames = c("all", "none")
## create an empty string of the input
inputString <<- ""
## specify the prompt message
wfoPromptMessage = paste("Which lists do you want to run a synonym check on? The possible commands are: \n lists names: ", paste(simpleListNames, collapse = ", "), "\n multi-action: ", paste(multiListNames, collapse = ", "), "\n")
## Run the command line input check
inputString = checkInputCommands(inputString, simpleListNames, multiListNames, wfoPromptMessage)

## Split input into individual commands
inputCommands = strsplit(inputString, ",")[[1]]
## Source the ABA list
if ("aba" %in% inputCommands || "all" %in% inputCommands) source("components/aba_wrangling.R")
if (!"aba" %in% inputCommands || "none" %in% inputCommands) cat("Skipping the sourcing of the ABA list")
## Source the AMBIO list
if ("ambio" %in% inputCommands || "all" %in% inputCommands) source("components/ambio_wrangling.R")
if (!"ambio" %in% inputCommands || "none" %in% inputCommands) cat("Skipping the sourcing of the AMBIO list")
## Source the regions and GBIF list
if ("gbif" %in% inputCommands || "all" %in% inputCommands) {
  source("components/regions_data_import.R")
  source("components/gbif_species_retrieval.R")
}
if (!"gbif" %in% inputCommands || "none" %in% inputCommands) cat("Skipping the sourcing of the GBIF list and region import")
## Source the GloNAF list
if ("glonaf" %in% inputCommands || "all" %in% inputCommands) source("components/glonaf_wrangling.R")
if (!"glonaf" %in% inputCommands || "none" %in% inputCommands) cat("Skipping the sourcing of the GloNAF list")

endTime = Sys.time()
### Calculate the time
formatted_elapsed_time = format_elapsed_time(starTime, endTime)
## Print message with elapsed time
cat("Sourcing scripts finished in ", formatted_elapsed_time, "\n")

## Source the synonym check
source("components/synonym_check.R")
cat("sourcing completed. \n")

## If synonym_check is not conducted,load files
if (!"aba" %in% inputCommands || "none" %in% inputCommands) {
  ### Load ABA
  wfo_aba_arctic_present = read.csv("outputs/wfo_aba_arctic_present.csv")
  wfo_aba_arctic_absent = read.csv("outputs/wfo_aba_arctic_absent.csv")
}
if (!"ambio" %in% inputCommands || "none" %in% inputCommands) {
  ### Load AMBIO
  wfo_ambio_arctic_present = read.csv("outputs/wfo_ambio_arctic_present.csv")
  wfo_ambio_arctic_absent = read.csv("outputs/wfo_ambio_arctic_absent.csv")
}
if (!"gbif" %in% inputCommands || "none" %in% inputCommands) {
  ### Load GBIF
  wfo_gbif_species = read.csv("outputs/wfo_gbif_species.csv")
}
if (!"glonaf" %in% inputCommands || "none" %in% inputCommands) {
  ### Load GloNAF
  wfo_glonaf_species = read.csv("outputs/wfo_glonaf_species.csv")
}

# choose Scientifically accepted names
message("---------- Selecting and removing duplicates from the ABA list ----------")
## ABA Arctic present
## Use WFO.one to remove synonyms of the same species to only be left with one name of the same species
aba_arctic_present = WFO.one(wfo_aba_arctic_present)
### Select the scientific names and remove NA
aba_arctic_present = aba_arctic_present %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All aba_arctic_present NAs removed:", any(!is.na(aba_arctic_present)))
## Remove duplicates
aba_arctic_present = distinct(aba_arctic_present)
cat("All aba_arctic_present are distinct:", any(!duplicated(aba_arctic_present)))

## ABA Arctic absent
## Use WFO.one to remove synonyms of the same species to only be left with one name of the same species
aba_arctic_absent = WFO.one(wfo_aba_arctic_absent)
##Select the scientific names
aba_arctic_absent = aba_arctic_absent %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All aba_arctic_absent NAs removed:", any(!is.na(aba_arctic_absent)))
## Remove duplicates
aba_arctic_absent = distinct(aba_arctic_absent)
cat("All aba_arctic_absent are distinct:", any(!duplicated(aba_arctic_absent)))

message("---------- Selecting and removing duplicates from the AMBIO list ----------")
## AMBIO Arctic present
## Use WFO.one to remove synonyms of the same species to only be left with one name of the same species
ambio_arctic_present = WFO.one(wfo_ambio_arctic_present)
### Select the scientific names and remove NA
ambio_arctic_present = ambio_arctic_present %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All ambio_arctic_present NAs removed:", any(!is.na(ambio_arctic_present)))
## Remove duplicates
ambio_arctic_present = distinct(ambio_arctic_present)
cat("All ambio_arctic_present are distinct:", any(!duplicated(ambio_arctic_present)))

## ambio Arctic absent
ambio_arctic_absent = WFO.one(wfo_ambio_arctic_absent)
ambio_arctic_absent = ambio_arctic_absent %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All ambio_arctic_absent NAs removed:", any(!is.na(ambio_arctic_absent)))
## Remove duplicates
ambio_arctic_absent = distinct(ambio_arctic_absent)
cat("All ambio_arctic_absent are distinct:", any(!duplicated(ambio_arctic_absent)))

# Combine AMBIO and ABA lists
message("------ Merging ABA and AMBIO lists ------")
## merge present lists
arctic_present = aba_arctic_present %>% 
  full_join(ambio_arctic_present, by = "scientificName") %>% 
  distinct(scientificName, .keep_all = TRUE)

## Run NA and distinct check
cat("All arctic_present NAs removed:", any(!is.na(arctic_present)))
cat("All arctic_present values are distinct: ", any(!duplicated(arctic_present)))

## merge absent lists
arctic_absent = aba_arctic_absent %>% 
  full_join(ambio_arctic_absent, by = "scientificName") %>% 
  distinct(scientificName, .keep_all = TRUE)

## Run NA and distinct check
cat("All arctic_absent NAs removed:", any(!is.na(arctic_absent)))
cat("All arctic_absent values are distinct: ", any(!duplicated(arctic_absent)))
cat("Merging of ABA and AMBIO finished")

message("------ Selecting and removing duplicates from the GBIF list ------")

## GBIF species
gbif_species = WFO.one(wfo_gbif_species)
gbif_species = gbif_species %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All gbif_species NAs removed:", any(!is.na(gbif_species)))
## Remove duplicates
gbif_species = distinct(gbif_species)
cat("All gbif_species are distinct:", any(!duplicated(gbif_species)))

# Merge GBIF list with ABA and AMBIO list
message("------ Merging the combined ABA and AMBIO list with the GBIF list ------")
## Remove rows from filtered_species data frame that have matching values to the arctic_absent species. This is to avoid duplicates.
gbif_arctic_present = merge(gbif_species, 
                                arctic_present, 
                                by = "scientificName"
)
## Remove rows from gbif_species data frame that have matching values to the arctic_present species
gbif_arctic_absent = anti_join(gbif_species, 
                         arctic_present, 
                         by = "scientificName"
)

## Remove rows that are similar to the arctic_absent to avoid duplicates
gbif_arctic_absent = anti_join(gbif_arctic_absent, 
                               arctic_absent, 
                               by = "scientificName"
)
## Check how many are similar
cat("Number of rows in GBIF that are similar to the arctic_absent list: ", nrow(anti_join(gbif_species, arctic_present, by = "scientificName")) - nrow(anti_join(gbif_arctic_absent, arctic_absent, by = "scientificName")), "\n")

## Now add all the absent species to the list. 
gbif_arctic_absent = bind_rows(gbif_arctic_absent, arctic_absent)

cat("gbif_arctic_present and gbif_arctic_absent lists have been made")

# Merge filtered species with the GloNAF list
message("------ Selecting and removing duplicates from the GloNAF list ------")
## Glonaf species
glonaf_species = WFO.one(wfo_glonaf_species)
glonaf_species = glonaf_species %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))

## Remove rows from glonaf_species data frame that have matching values to the arctic_present species
glonaf_arctic_present = merge(glonaf_species, 
                                         arctic_present, 
                                         by = "scientificName"
)

## Remove rows from glonaf_species data frame that have matching values to the arctic_present species
## This one will merge with arctic_absent species later
glonaf_arctic_absent = anti_join(glonaf_species, 
                                         arctic_present, 
                                         by = "scientificName"
)

cat("All glonaf_arctic_absent NAs removed:", any(!is.na(glonaf_arctic_absent)))
## Remove duplicates
glonaf_arctic_absent = distinct(glonaf_arctic_absent)
cat("All glonaf_arctic_absent are distinct:", any(!duplicated(glonaf_arctic_absent)))

# Merge the filtered list with the GloNAF list
message("------ Merging the gbif_arctic_absent list with the glonaf_arctic_absent list ------")
## Remove the glonaf species from the filtered_list. This is to avoid duplicates when merging
filtered_species = anti_join(glonaf_arctic_absent, 
                             gbif_arctic_absent,
                             by = "scientificName"
                             )
## Check how many are similar
cat("Number of rows in gbif_arctic_absent that are similar to the glonaf_arctic_absent list: ", nrow(merge(glonaf_arctic_absent, gbif_arctic_absent)), "\n")

## Merge filtered_species from glonaf with the gbif_species
filtered_species = bind_rows(gbif_arctic_absent, filtered_species)
## Check for NAs and duplications
cat("All filtered_species NAs removed:", any(!is.na(filtered_species)))
cat("All filtered_species are distinct:", any(!duplicated(filtered_species)))
cat("glonaf_arctic_absent has been merged with gbif_arctic_absent")

## Replace odd multiplication sign with x in order to better display the final list
filtered_species$scientificName = lapply(filtered_species$scientificName, function(x) gsub("(\\w) × (\\w)", "\\1 x \\2", x))
filtered_species$scientificName = unlist(filtered_species$scientificName)

message("----- Running similarity check -----")
## Run a similarity check to find possible missed duplications
similarityCheck_filtered_species = similarity_check(filtered_species, "scientificName", "scientificName", "jw", 0.05)
similarityCheck_filtered_species = data.frame(similarityCheck_filtered_species)

## Filter out different variances
similarityCheck_filtered_species = filter_rows_after_split_text(similarityCheck_filtered_species, "name", "similarName", "var.")
## Filter out hybrids
similarityCheck_filtered_species = filter_rows_around_split_text(similarityCheck_filtered_species, "name", "similarName", "x")
## Filter out different subspecies
similarityCheck_filtered_species = filter_rows_after_split_text(similarityCheck_filtered_species, "name", "similarName", "subsp.")

## Use rgbif to check the similar names
## combine the two columns into one
scientific_names = unlist(similarityCheck_filtered_species[ ,c("name", "similarName")] )

## Only look in the vascular plants
higherTaxonKey = name_backbone(name = "Tracheophyta")$usageKey

speciesKeys = lapply(scientific_names, function(x) {
  name_lookup(x, higherTaxonKey = higherTaxonKey)$data$speciesKey
})

## Find the maximum number of speciesKey values
max_speciesKeys = max(sapply(speciesKeys, function(x) length(unique(x))))
## Create matrix
mat = matrix(nrow = length(speciesKeys), ncol = max_speciesKeys + 1)
mat[, 1] = scientific_names
for (i in seq_along(speciesKeys)) {
  mat[i, 2:(length(unique(speciesKeys[[i]])) + 1)] = unique(speciesKeys[[i]])
}
## Convert matrix to data frame
speciesKeys_unique = as.data.frame(mat, stringAsFactors = FALSE)
## Convert the value columns to numeric
speciesKeys_unique[, -1] = lapply(speciesKeys_unique[, -1], as.numeric)
## Set the column names
colnames(speciesKeys_unique) = c("scientificName", "speciesKey")
## remove all other columns to only keep one value
speciesKeys_unique = speciesKeys_unique[, c(1, 2)]
## Write it out
write.csv(speciesKeys_unique, "outputs/similarityCheck_filtered_species.csv", fileEncoding = "UTF-8")
## Get duplicated species
duplicate_species = speciesKeys_unique[duplicated(speciesKeys_unique$speciesKey), ]
## Filter duplicated species
filtered_duplicate_species = speciesKeys_unique[!duplicated(speciesKeys_unique$speciesKey), ]

cat("Filtering out", nrow(duplicate_species), "Species from the filtered_species list: \n", paste(duplicate_species$scientificName, collapse = "\n"), "\n")

## Make the final list
filtered_species_final = anti_join(filtered_species, duplicate_species, by = "scientificName")

cat("------ Writing the final filtered lists to CSV ------")
## Write out the final CSV of all species outside the Arctic
write.csv(filtered_species_final, "outputs/filtered_species_final.csv", fileEncoding = "UTF-8")
end_filterScript_time = Sys.time()
### Calculate the time
formatted_elapsed_time = format_elapsed_time(start_filterScript_time, end_filterScript_time)
## Print message with elapsed time
message("Sourcing scripts finished in ", formatted_elapsed_time, "\n")
