source("./src/utils/utils.R")
source("./src/setup/setup.R")
source("./src/filter/filter.R")
source("./src/hypervolume/hypervolume.R")

sp_occ <- filter_process(
  test = NULL,
  cores.max = 14,
  verbose = T
  )

#sp_df <- setup_sp(test = "small")

min_disk_space <- get_disk_space("/home", units = "GB") * 0.2

## Potential issue where it is not seperating the species properly and running all of them on every core. FIX

# Check memory peak of one node by conducting a test run
peak_ram <- peakRAM({
  parallell_processing(
    sp_list, 
    method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
    accuracy = "accurate",
    project = "laea",
    proj.incl.t = 0.5,
    iterations = NULL,
    max_cores = 5,
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
  method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
  accuracy = "fast",
  project = "laea",
  proj.incl.t = c(0.5, 0.25, 0.1),
  iterations = NULL,
  max_cores = max_cores, 
  min.disk.space = min_disk_space,
  show.plot = F,
  verbose = F
)


t <- fread("./outputs/filter/arctic/chunk/species/A")

t1 <- fread("./outputs/filter/arctic/sp_w_keys.csv")

any(is.na(t1))
length(which(is.na(t1)))

t3 <- fread("./outputs/filter/test/test-big/chunk/species/Calypso-bulbosa.csv")
View(t3)
