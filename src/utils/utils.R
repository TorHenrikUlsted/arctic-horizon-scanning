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
  "ggnewscale",
  "ggrepel",
  "ggpubr",
  "gridExtra",
  "extrafont",
  "tidyterra",
  "scales",
  "rnaturalearth",
  "rnaturalearthdata",
  "countrycode",
  "ggalluvial",
  "heatmaply",
  "plotrix",
  "viridis",
  "stringdist",
  "stringr",
  "knitr",
  "httr",
  "jsonlite"
)

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (Sys.info()["sysname"] == "Linux") {
  system.running <- "linux"
} else if (Sys.info()["sysname"] == "Darwin") {
  system.running <- "mac"
} else if (Sys.info()["sysname"] == "Windows") {
  system.running <- "windows"
} else {
  stop("Cannot determine what system is being used, add manually. e.g: system.running <- 'Linux'")
}

source("./src/utils/components/check_updates.R")
updated <- check_updates(pkgs)

if (updated) {
  warning("Session needs to be restarted due to package updates or installations.\n")
  
  # Check if rstudioapi is available
  if ("package:rstudioapi" %in% search()) {
    rm(ls())
    library(rstudioapi)
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

cat("Loading config handler. \n")
source("./src/utils/components/config_handler.R")

cat("Loading WGSRPD. \n")
download_github_dir_if("tdwg", "wgsrpd", "master", "level2", "resources/region/wgsrpd/level2")
download_github_dir_if("tdwg", "wgsrpd", "master", "level3", "resources/region/wgsrpd/level3")
