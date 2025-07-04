load_utils <- function(parallel = FALSE) {
  tryCatch({
    setwd("arctic-horizon-scanning")
  }, error = function(e) {
    cat("Working directory: ", getwd(), "\n")
  })
  
  pkgs = list(
    common = c(
      "tools",
      "data.table",
      "yaml",
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
      "patchwork",
      "gamlss",
      "mgcv",
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
      "vegan",
      "quantreg",
      "splines"
    ),
    non_parallel = c(
      "cli",
      "rgbif",
      "httr",
      "jsonlite",
      "rcrossref"
    )
  )
  
  if (!parallel) {
    pkgs <- c(pkgs$common, pkgs$non_parallel)
  } else {
    pkgs <- pkgs$common
  }
  
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  
  sys.source("./R/utils/components/check_updates.R", envir = globalenv())
  updated <- check_updates(pkgs)
  
  for (pkg in pkgs) {
    if (!pkg %in% (.packages())) {
      library(pkg, character.only = TRUE)
      cat("Package", pkg, "loaded into global environment.\n")
    }
  }
  
  if (updated) {
    warning("Session needs to be restarted due to package updates or installations.\n")
    
    # Check if rstudioapi is available
    if (!parallel && "package:rstudioapi" %in% search()) {
      rm(ls())
      library(rstudioapi)
      rstudioapi::restartSession()
    } else {
      message("Please restart your R session manually.\n")
    }
  }
  
  sys.source("./R/utils/components/custom_colors.R", envir = globalenv())
  assign("cc", custom_colors(), envir = globalenv())
  
  component_files <- c(
    
    "./R/utils/components/handlers/condition_handlers.R",
    "./R/utils/components/handlers/memory_handlers.R",
    "./R/utils/components/helpers/citation_helpers.R",
    "./R/utils/components/helpers/data_table_helpers.R",
    "./R/utils/components/helpers/ggplot_helpers.R",
    "./R/utils/components/helpers/object_helpers.R",
    "./R/utils/components/helpers/spatial_helpers.R",
    "./R/utils/components/helpers/system_helpers.R",
    "./R/utils/components/handlers/time_handlers.R",
    "./R/utils/components/handlers/file_handlers.R",
    "./R/utils/components/handlers/lock_handlers.R",
    "./R/utils/components/handlers/markdown_handlers.R",
    "./R/utils/components/conversion.R",
    "./R/utils/components/loader.R",
    "./R/utils/components/config.R"
  )
  
  for (file in component_files) {
    cat("Sourcing", file, "\n")
    sys.source(file, envir = globalenv())
  }
  
  # Initialize variables in global environment
  assign("WFO_file", load_wfo(), envir = globalenv())
  
  cat("Loading WGSRPD. \n")
  download_github_dir_if("tdwg", "wgsrpd", "master", "level2", "resources/region/wgsrpd/level2")
  download_github_dir_if("tdwg", "wgsrpd", "master", "level3", "resources/region/wgsrpd/level3")
  
  if (!parallel) {
    sys.source("./R/utils/components/handlers/cite_handlers.R", envir = globalenv())
    format <- "bibtex"
    
    if(!if_file(paste0("./outputs/utils/citations/citations.", format), "today")) {
      cat("Creating reference list. \n")
      cited_packages <- cite_packages(pkgs, formats = format)
      
      cat("Citing loaders. \n")
      cited_packages[[format]] <- c(cited_packages[[format]], load_apg(cite = TRUE))
      cited_packages[[format]] <- c(cited_packages[[format]], load_ppg(cite = TRUE))
      cited_packages[[format]] <- c(cited_packages[[format]], load_gpg(cite = TRUE))
      
      # Functions relying on packages
      write_citations(cited_packages, "./outputs/utils/citations")
    } else {
      cat("Reference list already created today. \n")
    }
  }
}
