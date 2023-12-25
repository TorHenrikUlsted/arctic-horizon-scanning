source("./src/utils/utils.R")
source("./src/setup/setup.R")
source("./src/hypervolume/hypervolume.R")

sp_df <- setup_sp(test = T, big_test = F)

sp_data <- sp_df %>% 
  select(species, decimalLongitude, decimalLatitude, coordinateUncertaintyInMeters, coordinatePrecision, countryCode, stateProvince, year)

sp_data <- sp_data[order(sp_data$species), ]

sp_list <- split(sp_data, sp_data$species)

min_disk_space <- get_disk_space("/home", units = "GB") * 0.2

# Check memory peak of one node by conducting a test run
peak_ram <- peakRAM({
  parallell_processing(
    sp_list, 
    method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
    projection = "longlat",
    proj.incl.t = c(0.5, 0.25, 0.1),
    iterations = 2,
    max_cores = 1, 
    min.disk.space = min_disk_space,
    show.plot = F,
    verbose = F
  ) 
})
 
invisible(gc())

peak_ram_gb <- peak_ram[[4]] / 1024

peak_ram_gb_man <- 18

max_cores <- floor(mem_limit/1024^3 / peak_ram_gb_man)

cores_to_use <- min(length(sp_list), detectCores() / 2, max_cores)

parallell_processing(
  sp_list, 
  method = "gaussian", #box approx 13 min, gaussian 1 hours 10 minutes
  projection = "longlat",
  proj.incl.t = c(0.5, 0.25, 0.1),
  iterations = NULL,
  max_cores = max_cores, 
  min.disk.space = min_disk_space,
  show.plot = F,
  verbose = F
)
