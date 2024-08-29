## Example script
source("./src/utils/utils.R")
load_utils()
source("./src/setup/setup_sequence.R")
source("./src/filter/filter_sequence.R")
source("./src/hypervolume/parallel_hypervolume.R")
source("./src/hypervolume/node_hypervolume.R")
source("./src/hypervolume/hypervolume.R")
source("./src/visualize/visualize.R")
source("./src/main.R")

# Can run tests with
# spec.known = "test_known"
# spec.unknown = "test_small" OR "test_big"

main(
  spec.known = "test_known",
  spec.unknown = "test_big",
  coord.uncertainty = NULL,
  gbif.occ.region = NULL, # If wanting to download files within a shapefile -- converted to WKT
  download.key = "0186013-240321170329656",
  download.doi = "https://doi.org/10.15468/dl.awqjxw",
  hv.iterations = NULL,
  hv.method = "box",
  hv.accuracy = "accurate",
  hv.dims = c(18, 10, 3, 4),
  hv.incl.threshold = 0.5,
  vis.shape = "./outputs/setup/region/cavm-noice/cavm-noice.shp", # change to simply name
  vis.title = TRUE,
  vis.region.name = "the Arctic", 
  vis.subregion.name = "Floristic Province", 
  vis.composition.taxon = "order",
  vis.save.device = "svg",
  vis.save.unit = "px",
  plot.show = FALSE,
  verbose = TRUE,
  force.seq = NULL
)

create_derived_dataset(
  occurrences.dir = paste0(
    "./outputs/filter/", 
    "glonaf",
    "/chunk/species"
  ),
  verbose = FALSE
)

pack_repository(
  filename = "Horizon-Scanning-Repository",
  which.sequence = "all"
)