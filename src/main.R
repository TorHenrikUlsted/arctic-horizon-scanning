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
peak_ram <- setup_hv_process(min_disk_space)

# Get 60% of total cores and get a ratio of high and low memory core usage
total_cores <- detectCores() * 0.6
cores_high <- round(total_cores * (5 / 8))
cores_low <- total_cores - cores_high

mem_high_gb <- peak_ram$high #GB
mem_low_gb <- peak_ram$low #GB
total_mem_gb <- mem_high_gb + mem_low_gb

mem_core_high <- mem_high_gb / total_cores
mem_core_low <- mem_low_gb / total_cores

cores_max_total <- min(length(sp_list), floor((mem_limit / 1024^3) / total_cores))
cores_max_high <- min(length(sp_list), floor((mem_limit / 1024^3) * (5 / 8) / mem_high_gb), cores_high)
cores_max_low <- min(floor((mem_limit / 1024^3) * (3 / 8) / mem_low_gb), cores_max_total - cores_max_high)

parallel_processing(
  spec.list = sp_list, # list of strings
  method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
  accuracy = "accurate",
  hv.dims = c(18, 10, 3, 4),
  hv.projection = "laea",
  proj.incl.t = 0.5,
  iterations = NULL,
  cores.max.high = cores_max_high, 
  cores.max = cores_max_total,
  min.disk.space = min_disk_space,
  hv.dir = "./outputs/hypervolume/sequence",
  show.plot = F,
  verbose = T
) 

#as.numeric(gsub("node", "", readLines("outputs/hypervolume/sequence/logs/node-iterations.txt")))
