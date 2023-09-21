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

cat("Changing CAVM projections \n")
## Change projection to same as the bio variables
cavm = terra::project(cavm, "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +no_defs +units=m")
crs(cavm, proj = T, describe = T)


cat("Creating Boreal Ecoregion shape  \n")
## read shapefile of ecoregions
wwfecoRegions = vect("resources/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp")
## check projection
crs(wwfecoRegions, proj = T, describe = T)
## check summary
wwfecoRegions
## filter out all others but the boreal biome
boreal_biome = wwfecoRegions[wwfecoRegions$BIOME == 6, ]

# Combine WKT to download GBIF data from all regions
cat("Combining WKTs \n")
# Define the two extents
ext_cavm = ext(cavm)
ext_boreal = ext(boreal_biome)
# Find the min and max y values
min_y = min(ext_boreal[3], ext_cavm[3])
max_y = max(ext_boreal[4], ext_cavm[4])
# Define the max x values
min_x = -180
max_x = 180

# Create the combined polygon
combined_WKT = sprintf("POLYGON((%f %f, %f %f, %f %f, %f %f, %f %f))",  min_x, min_y, max_x, min_y, max_x, max_y, min_x, max_y, min_x, min_y)
# Create the combined polygon
combined_WKT = sprintf("POLYGON((%s %s,%s %s,%s %s,%s %s,%s %s))",  
                       formatC(min_x, format = "f", digits = 5), formatC(min_y, format = "f", digits = 5), 
                       formatC(max_x, format = "f", digits = 5), formatC(min_y, format = "f", digits = 5), 
                       formatC(max_x, format = "f", digits = 5), formatC(max_y, format = "f", digits = 5), 
                       formatC(min_x, format = "f", digits = 5), formatC(max_y, format = "f", digits = 5), 
                       formatC(min_x, format = "f", digits = 5), formatC(min_y, format = "f", digits = 5))
# Split the string into individual numbers
numbers = strsplit(combined_WKT, " ")[[1]]
# Remove trailing zeros from each number
numbers = sapply(numbers, function(x) {
  x = gsub("0+$", "", x)
  gsub("\\.$", "", x)
})
# Combine the numbers back into a single string
combined_WKT = paste(numbers, collapse=" ")
# Remove spaces after commas
combined_WKT = gsub(", ", ",", combined_WKT)
combined_WKT

write(combined_WKT, "outputs/gbif_retrieval_process/combined_WKT.txt")
cat("CAVM and Boreal ecoregions imported \n")