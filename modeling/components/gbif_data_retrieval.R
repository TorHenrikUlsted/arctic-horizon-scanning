library(rgbif)

# 1. Download list of unique species names in the Arctic region

## ----------------------------------------------- Arctic Data --------------------------------------------------------------

## Use the WKT polygon created from the CAVM region
cavmWKT = "POLYGON((-169.50929 55.79623,172.06954 55.79623,172.06954 83.62742,-169.50929 83.62742,-169.50929 55.79623))"

## Download data from GBIF for the specified area using occ_search, this is only for testing purposes
arcticUniqueSpecies = unique(occ_search(geometry = cavmWKT, limit = 1000)$data$species)


## This is for the real download event
## Prepare GBIF keys
taxonName = "Tracheophyta"
taxonKey = name_backbone(taxonName)$usageKey
## Test data download
occ_download_prep(pred("taxonKey", taxonKey), 
                  pred("hasGeospatialIssue", FALSE),
                  pred("hasCoordinate", TRUE),
                  pred("geometry", cavmWKT),
                  pred("occurrenceStatus","PRESENT"), #could there be a point to also include absent? then use some filterin here as well.
                  pred_notin("basisOfRecord", c("LIVING_SPECIMEN", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN")),
                  format = "SIMPLE_CSV"
)

## Download the data 
vascPlants = occ_download(pred("taxonKey", taxonKey), 
                          pred("hasGeospatialIssue", FALSE),
                          pred("hasCoordinate", TRUE),
                          pred("geometry", cavmWKT),
                          pred("occurrenceStatus","PRESENT"),
                          pred_notin("basisOfRecord", c("LIVING_SPECIMEN", "PRESERVED_SPECIMEN", "FOSSIL_SPECIMEN")),
                          format = "SIMPLE_CSV"
)
## check status
occ_download_wait(vascPlants)

## get the download Data and import to create dataframe
vascPlants_file = occ_download_get(vascPlants, path = "resources", overwrite = T)
vascPlants_df = occ_download_import(vascPlants_file, path = "resources" )

## USE THIS IF ALREADY DOWNLOADED
vascPlants_df = occ_download_import(as.download("resources/0083144-230224095556074.zip"))


## test to see which branch includes NA
gbif_test_na_branch = vascPlants_df %>% 
  group_by(kingdom, phylum, class, order, family, genus) %>% 
  summarise(species = species)

apply(gbif_test_na_branch, 2, function(x) 
  any(is.na(x)) | any(is.infinite(x)) 
)

## check if any lat or long are NA
with(vascPlants_df, any(is.na(decimalLatitude)))
with(vascPlants_df, any(is.na(decimalLongitude)))
with(vascPlants_df, any(decimalLongitude < -169.50929 | decimalLongitude > 172.06954))
with(vascPlants_df, any(decimalLatitude < 55.79623 | decimalLatitude > 83.62742))

#Make into spatial points
vascPlantsLongLat = cbind(vascPlants_df$decimalLongitude, vascPlants_df$decimalLatitude)
sp_occ = vect(vascPlantsLongLat, vascPlants_df, type="points", crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +no_defs +units=m")
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

# 2. Download list of unique species names in the Boreal region

## ----------------------------------------------- Boreal Data --------------------------------------------------------------



