tryCatch({
  setwd("arctic-horizon-scanning")
}, error = function(e) {
  cat("Working directory: ", getwd(), "\n")
})  


pkgs = c(
  "rgbif",
  "dplyr",
  "WorldFlora",
  "geodata",
  "terra",
  "tidyterra",
  "sf",
  "sp",
  "data.table",
  "CoordinateCleaner",
  "spThin",
  "reticulate",
  "stringdist",
  "stringr",
  "parallel",
  "progressr",
  "crayon",
  "corrplot",
  "hypervolume",
  "dynRB",
  "VennDiagram",
  "ggplot2",
  "alphahull",
  "rgl"
)

source("./src/utils/components/check_updates.R")
updated <- check_updates(pkgs)

if (updated) {
  cat(yellow("Session needs to be restarted due to package updates or installations.\n"))
  
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

cat("Include the longlat and cavm laea CRS \n")

laea_crs <- "+proj=laea +lon_0=180 +lat_0=90 +datum=WGS84"

longlat_crs <- "+proj=longlat +datum=WGS84 +ellps=WGS84"

cat("Loading WFO file. \n")
source("./src/utils/components/get_wfo_backbone.R")

cat("Loading WWF Ecoregion file. \n")
source("./src/utils/components/get_wwf_ecoregions.R")

cat("Loading time tracker. \n")
source("./src/utils/components/time_tracker.R")

cat("Loading utf8 function. \n")
source("./src/utils/components/set_df_utf8.R")

cat("Loading duplicate logger. \n")
source("./src/utils/components/log_duplicates.R")

cat("Loading source all function. \n")
source("./src/utils/components/source_all.R")

cat("Loading similarity check function. \n")
source("./src/utils/components/stringdist_similarity_check.R")

cat("Loading create dir if \n")
source("./src/utils/components/create_dir_if.R")

cat("Loading input command check. \n")
source("./src/utils/components/check_input_cmd.R")

cat("Loading filter rows after split text. \n")
source("./src/utils/components/filter_rows_after_split_txt.R")

cat("Loading filter rows around split text. \n")
source("./src/utils/components/filter_rows_around_split_txt.R")

cat("Loading extract name after prefix. \n")
source("./src/utils/components/extract_name_after_prefix.R")
