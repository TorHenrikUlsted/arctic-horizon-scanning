## Example script
source("./src/utils/utils.R")
load_utils()
config$simulation$example = TRUE
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


# Need to edit the number of species, not 16
# Need to make sure the scientificName, is the accepted species and not synonym in row shorteing
# edit Error in if (is.na(highest_it)) 1 else highest_it + 1 : argument is of length zero
# We want to pass 1 species at a time, not a chunk. We only want the chunk to be made, then uncchunked in the same chunked order

main(
  spec.known = "arctic", # Name of combined present data 
  spec.unknown = "glonaf", # Name of combined absent data
  coord.uncertainty = NULL, # Can specify, else tries to estimate output resolution
  gbif.occ.region = NULL, # If wanting to download files within a shapefile -- converted to WKT
  download.key = "0033220-240906103802322",
  download.doi = "https://doi.org/10.15468/dl.fs9syg",
  hv.iterations = NULL,
  hv.method = "box",
  hv.accuracy = "accurate",
  hv.dims = c(18, 10, 3, 4),
  hv.incl.threshold = 0.5,
  vis.shape = "cavm-noice",
  vis.title = FALSE,
  vis.region.name = "the Arctic", # The name that will be displayed on plot titles
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
    "glonaf",
    "/chunk/species"
  ),
  verbose = FALSE
)

pack_repository(
  filename = "Horizon-Scanning-Repository",
  which.sequence = "all"
)
