## Example script
source("./src/utils/utils.R")
config$run$example = TRUE
source("./src/setup/setup_sequence.R")
source("./src/filter/filter_sequence.R")
source("./src/hypervolume/parallel_hypervolume.R")
source("./src/hypervolume/node_hypervolume.R")
source("./src/hypervolume/hypervolume.R")
source("./src/visualize/visualize.R")
source("./src/main.R")

main(
  spec.known = filter_arctic,
  spec.unknown = filter_glonaf,
  test = "small",
  approach = "precautionary",
  coord.uncertainty = NULL,
  climate.database = "worldclim",
  region = NULL,
  download.key = "0186013-240321170329656",
  download.doi = "https://doi.org/10.15468/dl.awqjxw",
  hv.iterations = NULL,
  hv.method = "box",
  hv.accuracy = "accurate",
  hv.dims = c(18, 10, 3, 4),
  hv.incl.threshold = 0.5,
  vis.shape = "./outputs/setup/region/cavm-noice/cavm-noice.shp",
  vis.projection = "laea",
  vis.title = TRUE,
  vis.region.name = "the Arctic", 
  vis.subregion.name = "Floristic Province", 
  vis.composition.taxon = "order",
  vis.save.device = "svg", 
  vis.save.unit = "px",
  plot.show = FALSE,
  verbose = TRUE,
  example = TRUE,
  force.seq = "setup"
)

create_derived_dataset(
  occurrences.dir = paste0(
    "./outputs/filter/", 
    gsub("filter_", "", get_obj_name(filter_glonaf)), 
    "/chunk/species"
  ),
  verbose = FALSE
)

pack_repository(
  filename = "Horizon-Scanning-Repository",
  which.sequence = "all"
)