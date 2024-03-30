source_all("./src/setup/components")
source_all("./src/setup/custom_setup")

setup_sequence <- function(cores.max = 1, verbose = FALSE) {
  
  wfo_speed <- check_system_speed(
    df.path = "./resources/data-raw/speed-test-species.csv",
    test.name = "wfo-speed",
    sample.size = NULL,
    cores.max = 1,
    fun = wfo_speed,
    verbose = verbose
  )
  
  system.speed.wfo <<- wfo_speed
  
  setup_raw_data(
    column = "rawName",
    test = NULL, 
    max.cores = cores.max,
    verbose = verbose,
    counter = 10
  )
  
  setup_region()
  
  setup_region_hv()
  
  setup_hv_sequence()
}