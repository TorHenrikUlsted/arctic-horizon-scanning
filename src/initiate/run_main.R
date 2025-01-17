source("./src/setup/setup_sequence.R")
source("./src/filter/filter_sequence.R")
source("./src/hypervolume/parallel_hypervolume.R")
source("./src/hypervolume/node_hypervolume.R")
source("./src/hypervolume/hypervolume.R")
source("./src/visualize/visualize.R")
source("./src/initiate/main.R")

# Can run tests with
# spec.known = "test_known"
# spec.unknown = "test_small" OR "test_big"

main(
  spec.known = config$dataset$known, # Name of combined present data
  spec.known.key = config$gbif$known$download.key,
  spec.known.doi = config$gbif$known$download.doi,
  spec.unknown = config$dataset$unknown, # Name of combined absent data
  spec.unknown.key = config$gbif$unknown$download.key,
  spec.unknown.doi = config$gbif$unknown$download.doi,
  gbif.occ.region = config$gbif$region, # If wanting to download files within a shapefile -- converted to WKT
  coord.uncertainty = config$projection$raster_scale_m,
  hv.iterations = config$hypervolume$iterations,
  hv.method = config$hypervolume$method,
  hv.accuracy = config$hypervolume$accuracy,
  hv.dims = config$hypervolume$dimensions,
  hv.incl.threshold = config$hypervolume$inclusion.threshold,
  vis.shape = config$visualization$shape,
  vis.title = config$visualization$title,
  vis.region.name = config$visualization$region.name, # The name that will be displayed on plot titles
  vis.projection = config$simulation$projection,
  vis.subregion.name = config$visualization$subregion.name, 
  vis.composition.taxon = config$visualization$composition.taxon, 
  vis.save.device = config$visualization$save.device,
  vis.save.unit = config$visualization$save.unit,
  plot.show = config$visualization$plot.show,
  validation = config$simulation$validation,
  total.cores = config$memory$total_cores,
  verbose = config$simulation$verbose,
  force.seq = config$simulation$force.sequence
)
