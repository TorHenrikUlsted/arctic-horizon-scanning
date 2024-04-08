source("./src/utils/utils.R")
source("./src/setup/setup_sequence.R")
source("./src/filter/filter_sequence.R")
source("./src/hypervolume/parallel_hypervolume.R")
source("./src/hypervolume/node_hypervolume.R")
source("./src/hypervolume/hypervolume.R")
source("./src/visualize/visualize.R")
source("./src/main.R")

main(
  spec.known = NULL,
  speck.unknown,
  test = NULL,
  column = "scientificName",
  coord.uncertainty = 4600,
  region = NULL,
  download.key = NULL,
  download.doi = NULL,
  hv.iterations = NULL,
  hv.method = "box",
  hv.accuracy = "accurate",
  hv.dims = NULL,
  hv.incl.threshold = 0.5,
  verbose = FALSE
)