message("Initiating hypervolume sequence")

source("utils.R")

shapefiles = c(
  cavm = "./resources/cavm2003/cavm.shp",
  wwfEcoRegion = "./resources/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp"
)
projection = crs("+proj=longlat +datum=WGS84")

source("./hypervolume/atoms/import_regions.R")
regions = import_regions(shapefiles, prj = projection, log_output = "./hypervolume/log.txt")
plot(regions$cavm)

# Get anticlockwise wkt (GBIF friendly)
source("./hypervolume/atoms/combine_wkt_anticlockwise.R")
anticlockwise_wkt = combine_wkt_anticlockwise(regions, max_x = T, min_x = T)

# Use long lat to get WorldClim data
source("./hypervolume/atoms/wc_crop_region.R")
wc_region = wc_to_region(regions$cavm)

source("./hypervolume/atoms/correlation.R")
biovars_importance <- get_correlation(wc_region)

# Pick the most important bioVariables
biovars_importance_sel <- biovars_importance[c(1, 3,4,5), ]