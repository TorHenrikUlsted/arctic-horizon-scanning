library(terra)

# Arctic CAVM
## Citation:
### CAVM Team. 2003. Circumpolar Arctic Vegetation Map. (1:7,500,000 scale), Conservation of Arctic Flora and Fauna (CAFF) 
### Map No. 1. U.S. Fish and Wildlife Service, Anchorage, Alaska. ISBN: 0-9767525-0-6, ISBN-13: 978-0-9767525-0-9

## Write SpatVector file of cavm shape
cavm = vect("modeling/resources/Cavm2003/aga_circumpolar_geobotanical_2003.shp")
crs(cavm, proj = T, describe = T)
plot(cavm)

## Change projection to same as the bio variables
cavm = terra::project(cavm, "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +no_defs +units=m")
crs(cavm, proj = T, describe = T)
plot(cavm)

# WWF Ecoregions --> Boreal
## read shapefile of ecoregions
wwfecoRegions = vect("modeling/resources/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp")
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
