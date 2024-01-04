source("./src/utils/utils.R")
source("./src/setup/setup.R")
source("./src/hypervolume/hypervolume.R")
shapefiles <- c(
  cavm = "./resources/region/cavm-noice/cavm-noice.shp"
)
regions <- import_regions(shapefiles, "./outputs/data_acquisition/region/logs")
setup_region_hv(regions$cavm, "cavm_noice", method = "box")
