sp_dir <- filter_process(
  test = NULL,
  cores.max = 12,
  verbose = T
  )

# Get file_names as a list
cat("Listing files in the directory:", cc$lightSteelBlue(sp_dir), "\n")

sp_list <- list.files(sp_dir, full.names = TRUE)

cat("Found", cc$lightSteelBlue(length(sp_list)), "species. \n")

min_disk_space <- get_disk_space("/export", units = "GB") * 0.2

# Check memory peak of one node by conducting a test run as well as setting up the entire hypervolume sequence
peak_ram <- setup_hv_sequence(min_disk_space)

mem_high_gb <- peak_ram$high #/ 1024 Removed GB transformation as it seems to be more correct that it is in GB already
mem_max_high <- floor(mem_limit/1024^3 / mem_high_gb * 0.9)
cores_max_high <- min(length(sp_list), floor(detectCores() * 0.6), mem_max_high)

mem_low_gb <- peak_ram$low #/ 1024 Removed GB transformation as it seems to be more correct that it is in GB already
mem_max_low <- floor(mem_max_high / mem_low_gb)
if (cores_max_high < mem_max_high) cores_max_low <- 0 else cores_max_low <- floor(detectCores() * 0.2 - cores_max_high)

cores_max <- cores_max_high + cores_max_low

parallell_processing(
  spec.list = sp_list, # list of strings
  method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
  accuracy = "accurate",
  hv.projection = "laea",
  proj.incl.t = 0.5,
  iterations = NULL,
  cores.max.high = cores_max_high, 
  cores.max = cores_max,
  min.disk.space = min_disk_space,
  hv.dir = "./outputs/hypervolume",
  show.plot = F,
  verbose = T
) 
