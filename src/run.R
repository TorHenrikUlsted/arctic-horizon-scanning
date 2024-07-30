source("./src/utils/utils.R")
source("./src/setup/setup_sequence.R")
source("./src/filter/filter_sequence.R")
source("./src/hypervolume/parallel_hypervolume.R")
source("./src/hypervolume/node_hypervolume.R")
source("./src/hypervolume/hypervolume.R")
source("./src/visualize/visualize.R")
source("./src/main.R")

main(
  spec.known = filter_known,
  spec.unknown = filter_unknown,
  test = "small",
  approach = "precautionary",
  coord.uncertainty = NULL,
  region = NULL,
  download.key = NULL,
  download.doi = NULL,
  hv.iterations = NULL,
  hv.method = "box",
  hv.accuracy = "accurate",
  hv.dims = NULL,
  hv.incl.threshold = 0.5,
  vis.shape = "./outputs/setup/region/region-name/region-name.shp",
  vis.projection = "laea",
  vis.title = TRUE,
  vis.region.name = "Region name", 
  vis.subregion.name = "Floristic Province", 
  vis.composition.taxon = "order",
  vis.save.device = "svg", 
  vis.save.unit = "px",
  plot.show = FALSE,
  verbose = FALSE,
  force.seq = NULL
)

create_derived_dataset(
  occurrences.dir = paste0(
    "./outputs/filter/", 
    gsub("filter_", "", get_obj_name(filter_unknown)), 
    "/chunk/species"
  ),
  verbose = FALSE
)

pack_repository(
  filename = "Horizon-Scanning-Repository",
  which.sequence = "all"
)
