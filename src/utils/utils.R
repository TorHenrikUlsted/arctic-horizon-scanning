tryCatch({
  setwd("arctic-horizon-scanning")
}, error = function(e) {
  cat("Working directory: ", getwd(), "\n")
})  
  
pkgs = c(
  "data.table",
  "rgbif",
  "dplyr",
  "WorldFlora",
  "geodata",
  "terra",
  "CoordinateCleaner",
  "spThin",
  "hypervolume",
  "crayon",
  "corrplot",
  "parallel",
  "sf",
  "sp",
  "ggplot2",
  "tidyterra",
  "scales",
  "rnaturalearth",
  "rnaturalearthdata",
  "ggalluvial",
  "heatmaply",
  "plotrix",
  "viridis",
  "stringdist",
  "stringr",
  "progressr"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

source("./src/utils/components/check_updates.R")
updated <- check_updates(pkgs)

if (updated) {
  warning("Session needs to be restarted due to package updates or installations.\n")
  
  # Check if rstudioapi is available
  if ("package:rstudioapi" %in% search()) {
    rstudioapi::restartSession()
  } else {
    message("Please restart your R session manually.\n")
  }
}

cat("Creating custom colors. \n")
source("./src/utils/components/custom_colors.R")
cc <- custom_colors()

cat("Creating reference list. \n")
source("./src/utils/components/cite_packages.R")
cite_packages(pkgs, formats = "bibtex")

cat("Loading condition handlers. \n")
source("./src/utils/components/condition_handlers.R")

cat("Loading memory handlers \n")
source("./src/utils/components/memory_handlers.R")

cat("Loading helper functions.\n")
source("./src/utils/components/helper_functions.R")

cat("Loading time tracker. \n")
source("./src/utils/components/time_tracker.R")

cat("Loading file_managers. \n")
source("./src/utils/components/file_managers.R")

cat("Loading lock_file. \n")
source("./src/utils/components/lock_handlers.R")

cat("Loading conversions. \n")
source("./src/utils/components/conversion.R")

cat("Setting up loaders. \n")
source("./src/utils/components/loader.R")
WFO_file <- load_wfo()

cat("Setting up config. \n")
source("./src/utils/components/config.R")
