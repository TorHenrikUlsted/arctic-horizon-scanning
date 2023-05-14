library(rgbif)
source("components/regions_data_import.R")

## Set parameteres
tracheophyta = name_backbone(name = "Tracheophyta", rank = "phylum")

## Search for unique species names within the tracheophyta phylum and the given geometry
gbif_unique_species = occ_search(geometry = combined_WKT, taxonKey = tracheophyta$usageKey, limit = 5000)
message("gbif_unique_species is done searching")

## Extract the unique_species and make them into a dataframe
gbif_unique_species = unique(data.frame(gbif_unique_species$data$species))

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
# gbif_species = occ_download(pred("taxonKey", taxonKey), 
#                                pred("hasGeospatialIssue", FALSE),
#                                pred("hasCoordinate", TRUE),
#                                pred("geometry", combined_WKT),
#                                pred("occurrenceStatus","PRESENT"),
#                                pred_in("basisOfRecord", c("OBSERVATION", "MACHINE_OBSERVATION", "HUMAN_OBSERVATION", "MATERIAL_SAMPLE", "OCCURENCE")),
#                                format = "SIMPLE_CSV"
# )

## check status
occ_download_wait(gbif_species)

## get the download Data and import to create dataframe
gbif_species_file = occ_download_get(gbif_species, path = "resources", overwrite = T)
gbif_species_df = occ_download_import(gbif_species, path = "resources" )

## USE THIS IF ALREADY DOWNLOADED
gbif_species_df = occ_download_import(as.download("resources/0083144-230224095556074.zip"))


## test to see which branch includes NA
gbif_test_na_branch = gbif_species_df %>% 
  group_by(kingdom, phylum, class, order, family, genus) %>% 
  summarise(species = species)

apply(gbif_test_na_branch, 2, function(x) 
  any(is.na(x)) | any(is.infinite(x)) 
)

## check if any lat or long are NA
with(gbif_species_df, any(is.na(decimalLatitude)))
with(gbif_species_df, any(is.na(decimalLongitude)))
## this has changed
with(gbif_species_df, any(decimalLongitude < -169.50929 | decimalLongitude > 172.06954))
with(gbif_species_df, any(decimalLatitude < 55.79623 | decimalLatitude > 83.62742))

#Make into spatial points
gbif_species_df_LongLat = cbind(gbif_species_df$decimalLongitude, gbif_species_df$decimalLatitude)
sp_occ = vect(gbif_species_df_LongLat, gbif_species_df, type="points", crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +no_defs +units=m")
sp_occ

#test plot
plot(cavm)
plot(sp_occ, add=T, col="red")

#Crop GBIF points to Arctic CAVM
#take time -> ETC 43 minutes
startTime = Sys.time()

cropGBIF = crop(sp_occ, cavm)
plot(cropGBIF)
sp_occMask = mask(cropGBIF, cavm)
plot(sp_occMask)

#Check endTime
endTime = Sys.time()
#print the time it took to complete the function
message("Cropping GBIF data to the Arctic complete. Elapsed time: ", endTime - startTime)

#plot cropped data
plot(cavm)
plot(sp_occMask, add=T, col="red")


# ----------- old code ------------ #

#Make corpped data into dataframe from 41 900 entries to 10 096 entries
sp_occCavm_df = as.data.frame(sp_occMask)

#uniqe number of classes
n_distinct(sp_occCavm_df$class)
#uniqe number of orders
n_distinct(sp_occCavm_df$order)
#uniqe number of families
n_distinct(sp_occCavm_df$family)
#uniqe number of genuses
n_distinct(sp_occCavm_df$genus)
#uniqe number of species
n_distinct(sp_occCavm_df$species)

#check for empty strings or NA
any(unique(sp_occCavm_df$species == ""))
any(is.na(unique(sp_occCavm_df$species)))

#List of species occurrences in the CAVM
cavmSpList = as.data.frame(unique(sp_occCavm_df$species))
cavmSpList = cavmSpList[!(cavmSpList$`unique(sp_occCavm_df$species)` == ""), ]
cavmSpList = as.data.frame(cavmSpList)
cavmSpList = `colnames<-`(cavmSpList, "species")

#Check for empty strings on new dataframe
any(cavmSpList == "")

#write species list to CSV
write.csv(cavmSpList, "outputs/Species in the CAVM.csv", row.names = F)