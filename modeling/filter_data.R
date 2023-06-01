# Filter all species lists

# Source scripts ETC 10 min without the synonym check
start_time = Sys.time()
cat("Starting to source scripts. \n")
## Source utility script
source("components/utils.R")
## Source the regions
source("components/regions_data_import.R")
## Source the ABA list
source("components/aba_wrangling.R")
## Source the AMBIO list
source("components/ambio_wrangling.R")
## Source the GBIF list
source("components/gbif_species_retrieval.R")
## Source the GloNAF list
source("components/glonaf_wrangling.R")
end_time = Sys.time()
### Calculate the time
formatted_elapsed_time = format_elapsed_time(start_time, end_time)
## Print message with elapsed time
cat("Sourcing scripts finished in ", formatted_elapsed_time, "\n")
## Source the synonym check
source("components/synonym_check.R")
cat("sourcing completed. \n")

## If synonym_check is not conducted,load files
if ("none" %in% inputCommands) {
  cat("Loading CSV files. \n")
  ### Load ABA
  wfo_aba_arctic_present = read.csv("outputs/wfo_aba_arctic_present.csv")
  wfo_aba_arctic_absent = read.csv("outputs/wfo_aba_arctic_absent.csv")
  ### Load AMBIO
  wfo_ambio_arctic_present = read.csv("outputs/wfo_ambio_arctic_present.csv")
  wfo_ambio_arctic_absent = read.csv("outputs/wfo_ambio_arctic_absent.csv")
  ### Load GBIF
  wfo_gbif_species = read.csv("outputs/wfo_gbif_species.csv")
  ### Load GloNAF
  wfo_glonaf_species = read.csv("outputs/wfo_glonaf_species.csv")
}
# choose Scientifically accepted names
cat("---------- Selecting and removing duplicates from the ABA list ----------")
## ABA Arctic present
### Select the scientific names and remove NA
aba_arctic_present = wfo_aba_arctic_present %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All aba_arctic_present NAs removed:", any(!is.na(aba_arctic_present)))
## Remove duplicates
aba_arctic_present = distinct(aba_arctic_present)
cat("All aba_arctic_present are distinct:", any(!duplicated(aba_arctic_present)))

## ABA Arctic absent
aba_arctic_absent = wfo_aba_arctic_absent %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All aba_arctic_absent NAs removed:", any(!is.na(aba_arctic_absent)))
## Remove duplicates
aba_arctic_absent = distinct(aba_arctic_absent)
cat("All aba_arctic_absent are distinct:", any(!duplicated(aba_arctic_absent)))

cat("---------- Selecting and removing duplicates from the AMBIO list ----------")
## AMBIO Arctic present
### Select the scientific names and remove NA
ambio_arctic_present = wfo_ambio_arctic_present %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All ambio_arctic_present NAs removed:", any(!is.na(ambio_arctic_present)))
## Remove duplicates
ambio_arctic_present = distinct(ambio_arctic_present)
cat("All ambio_arctic_present are distinct:", any(!duplicated(ambio_arctic_present)))

## ambio Arctic absent
ambio_arctic_absent = wfo_ambio_arctic_absent %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All ambio_arctic_absent NAs removed:", any(!is.na(ambio_arctic_absent)))
## Remove duplicates
ambio_arctic_absent = distinct(ambio_arctic_absent)
cat("All ambio_arctic_absent are distinct:", any(!duplicated(ambio_arctic_absent)))

# Combine AMBIO and ABA lists
cat("------ Merging ABA and AMBIO lists ------")
## merge present lists
arctic_present = aba_arctic_present %>% 
  full_join(ambio_arctic_present, by = "scientificName") %>% 
  distinct(scientificName, .keep_all = TRUE)

## merge absent lists
arctic_absent = aba_arctic_absent %>% 
  full_join(ambio_arctic_absent, by = "scientificName") %>% 
  distinct(scientificName, .keep_all = TRUE)
cat("Merging of ABA and AMBIO finished")

cat("------ Selecting and removing duplicates from the GBIF list ------")

## GBIF species 
gbif_species = wfo_gbif_species %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All gbif_species NAs removed:", any(!is.na(gbif_species)))
## Remove duplicates
gbif_species = distinct(gbif_species)
cat("All gbif_species are distinct:", any(!duplicated(gbif_species)))

# Merge GBIF list with ABA and AMBIO list
cat("------ Merging the combined ABA and AMBIO list with the GBIF list to make a filtered_species list ------")
## Remove rows from gbif_species data frame that have matching values to the arctic_present species
filtered_species = anti_join(gbif_species, 
                         arctic_present, 
                         by = c("scientificName" = "scientificName")
)

## Remove rows from filtered_species data frame that have matching values to the arctic_absent species. This is to avoid duplicates.
filtered_species = anti_join(filtered_species, 
                          arctic_absent, 
                          by = c("scientificName" = "scientificName")
)
## Now add all the absent species to the list. 
filtered_species = bind_rows(filtered_species, arctic_absent)

cat("filtered_species list has been made")

# Merge filtered species with the GloNAF list
cat("------ Selecting and removing duplicates from the GloNAF list ------")
## Glonaf species 
glonaf_species = wfo_glonaf_species %>% 
  select(scientificName) %>% 
  filter(!is.na(scientificName))
cat("All glonaf_species NAs removed:", any(!is.na(glonaf_species)))
## Remove duplicates
glonaf_species = distinct(glonaf_species)
cat("All glonaf_species are distinct:", any(!duplicated(glonaf_species)))

# Merge the filtered list with the GloNAF list
cat("------ Merging the filtered_species list with the GloNAF list ------")
## Remove the glonaf species from the filtered_list. This is to avoid duplicates when merging
filtered_species = anti_join(filtered_species, 
                             glonaf_species,
                             by = c("scientificName" = "scientificName")
                             )

## Merge filtered_species with the glonaf_species
filtered_species = bind_rows(filtered_species, glonaf_species)
## Check for NAs and duplications
cat("All filtered_species NAs removed:", any(!is.na(filtered_species)))
cat("All filtered_species are distinct:", any(!duplicated(filtered_species)))
cat("glonaf_species has been merged with filtered_species")

## Replace odd multiplication sign with x
filtered_species$scientificName = lapply(filtered_species$scientificName, function(x) gsub("(\\w) × (\\w)", "\\1 x \\2", x))
filtered_species$scientificName = unlist(filtered_species$scientificName)

cat("----- Running similarity check -----")
## Run a similarity check to find possible missed duplications
similarityCheck_filtered_species = similarity_check(filtered_species, "scientificName", "scientificName", "jw", 0.05)
write.csv(similarityCheck_filtered_species, "outputs/similarityCheck_filtered_species.csv", fileEncoding = "UTF-8")

cat("------ Writing the final filtered list to CSV ------")
## Write out the final CSV of all species outside the Arctic.
write.csv(filtered_species, "outputs/filtered_species.csv", fileEncoding = "UTF-8")
