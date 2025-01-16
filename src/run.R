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
  spec.known = NULL,
  spec.unknown = NULL,
  gbif.occ.region = NULL, # If wanting to download files within a shapefile -- converted to WKT
  hv.iterations = NULL,
  hv.method = "box",
  hv.accuracy = "accurate",
  hv.dims = NULL,
  hv.incl.threshold = 0.5,
  vis.shape = "shape-name", # Name of shapefile without .shp
  vis.title = FALSE,
  vis.region.name = "Region", # The name that will be displayed on plot titles
  vis.subregion.name = "Floristic Province", 
  vis.composition.taxon = "order", 
  vis.save.device = "svg",
  vis.save.unit = "px",
  plot.show = FALSE,
  verbose = FALSE,
  force.seq = NULL
)

create_derived_dataset(
  data.name = list(
    spec.known = "",
    spec.unknown = ""
  ),
  verbose = FALSE
)

pack_repository(
  filename = "Horizon-Scanning-Repository",
  which.sequence = "all"
)
