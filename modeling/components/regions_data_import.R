message("Initiating region import")
# Arctic CAVM
## Citation:
### CAVM Team. 2003. Circumpolar Arctic Vegetation Map. (1:7,500,000 scale), Conservation of Arctic Flora and Fauna (CAFF) 
### Map No. 1. U.S. Fish and Wildlife Service, Anchorage, Alaska. ISBN: 0-9767525-0-6, ISBN-13: 978-0-9767525-0-9

cat("Creating CAVM shape \n")
## Write SpatVector file of cavm shape
cavm = vect("resources/cavm2003/aga_circumpolar_geobotanical_2003.shp")
getwd()
crs(cavm, proj = T, describe = T)
plot(cavm)

cat("Changing CAVM projections \n")
## Change projection to same as the bio variables
cavm = terra::project(cavm, "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +no_defs +units=m")
crs(cavm, proj = T, describe = T)
plot(cavm)

cat("Creating and reversing CAVM WKT \n")
## Create WKT which is counter-clockwise as per GBIF standard
cavmWKT = sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))", ext(cavm)[1], ext(cavm)[3], ext(cavm)[2], ext(cavm)[3], 
                  ext(cavm)[2], ext(cavm)[4], ext(cavm)[1], ext(cavm)[4], ext(cavm)[1], ext(cavm)[3])
cavmWKT

cat("Creating Boreal Ecoregion shape  \n")
# WWF Ecoregions --> Boreal
## read shapefile of ecoregions
wwfecoRegions = vect("resources/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp")
## check projection
crs(wwfecoRegions, proj = T, describe = T)
## Check plot
plot(wwfecoRegions)
## check summary
wwfecoRegions

## filter out all others but the boreal biome
boreal_biome = wwfecoRegions[wwfecoRegions$BIOME == 6, ]
## Check plot
plot(boreal_biome)
## check summary
boreal_biome

cat("Creating and reversing Boreal WKT \n")
borealWKT = sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))", ext(boreal_biome)[1], ext(boreal_biome)[3], ext(boreal_biome)[2], ext(boreal_biome)[3], 
                    ext(boreal_biome)[2], ext(boreal_biome)[4], ext(boreal_biome)[1], ext(boreal_biome)[4], ext(boreal_biome)[1], ext(boreal_biome)[3])
borealWKT

# Combine WKT to download GBIF data from all regions
cat("Combining WKTs \n")
## Convert the WKT strings to sfc objects
cavm_sfc = st_as_sfc(cavmWKT)
boreal_sfc = st_as_sfc(borealWKT)
## Use st_union to combine the sfc objects into a single sfc object
combined_sfc = st_union(cavm_sfc, boreal_sfc)

#combined_sfc = c(cavm_sfc, boreal_sfc)
## Convert the combined sfc object back to a WKT string
combined_WKT = st_as_text(combined_sfc)

combined_WKT

write(combined_WKT, "outputs/gbif_retrieval_process/combined_WKT.txt")

cat("CAVM and Boreal ecoregions imported \n")