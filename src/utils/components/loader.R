load_wfo <- function() {
  if (!file.exists("./resources/wfo/classification.csv")) {
    if (!file.exists("./resources/wfo/WFO_Backbone.zip")) {
      tryCatch({
        if (!dir.exists("./resources/wfo")) dir.create("./resources/wfo", recursive = T)
        cat(blue("A progress window pops up here, check your taskbar \n"))
        suppressWarnings(WFO.download(save.dir = "./resources/wfo", WFO.remember = TRUE, timeout = 2000))
      }, error = function(e) {
        cat(red("Error in download:"), e$message, "\n")
        cat(red("Could not download WFO Backbone. Manual download needed. \n"))
        
        cat("Opening download page \n")
        browseURL("https://www.worldfloraonline.org/downloadData;jsessionid=D1501051E49AE20AB4B7297D021D6324")
      }, finally = {
        WFO_file <- "./resources/wfo/classification.csv"
      })
    } else {
      unzip("./resources/wfo/WFO_Backbone.zip", exdir = "./resources/wfo")
      WFO_file <- "./resources/wfo/classification.csv"
    }
  } else {
    WFO_file <- "./resources/wfo/classification.csv"
  }
  
  return(WFO_file)
}

load_wwf_ecoregions <- function() {
  if(!file.exists("resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp")) {
    downld_url = "https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip"
    fallbck_url = "https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world"
    paper_url = "https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2"
    zip_file = "resources/region/wwfTerrestrialEcoRegions.zip"
    dest_dir = "resources/region/wwfTerrestrialEcoRegions"
    official_dir = file.path("resources/region/wwfTerrestrialEcoRegions", "official")
    
    if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = T)
    
    message("WWF Terrestrial EcoRegions does not exist. trying to download...")
    
    tryCatch({
      ## Try to download the file from the primary URL, direct download
      download.file(downld_url, zip_file)
      ## Unzip the downloaded file
      unzip(zip_file, exdir = dest_dir)
      ## move the files to the desired directory
      files = list.files(official_dir, full.names = TRUE)
      file.rename(files, file.path(dest_dir, basename(files)))
      ## Clean up
      unlink(official_dir, recursive = TRUE)
      file.remove(zip_file)
    }, error = function(e) {
      message("Could not find the download file. trying to open website")
      tryCatch({
        ## Open the file in your default web browser
        browseURL(fallbck_url)
      }, error = function(e) {
        message("WebBrowser also failed, trying last resort. Opening the scientific paper.")
        browseURL(paper_url)
      })
    })
  } else {
    cat("WWF Ecoregions found. \n")
  }
}

load_region <- function(region.file, verbose = FALSE) {
  vebcat("Loading region", color = "funInit")
  
  if (!is.character(region.file)) {
    vebcat("Region is not a filepath.", color = "fatalError")
    stop()
  }
  
  prj = projection
  
  # Load the file using either rast or vect based on the file extension
  if (grepl("\\.shp$", region.file)) {
    region <- terra::vect(region.file)
  } else if (grepl("\\.tif$", region.file) || grepl("\\.nc$", region.file)) {
    region <- terra::rast(region.file)
  } else {
    stop("Unsupported file type", region.file)
  }
  
  # Get and print the extent
  ext_region <- ext(region)
  original_crs <- crs(region, proj = TRUE)
  
  ext_region_printable <- as.vector(ext_region)
  vebcat("Extent of region:", ext_region_printable, veb = verbose)
  
  # Check the original CRS
  vebcat("Original CRS: ", original_crs, veb = verbose)
  
  # If the original CRS is not correctly defined, define it
  if (is.na(original_crs) || original_crs == "") {
    vebcat("Found blank or na crs. Needs manual processing", color = "nonFatalError")
  }
  
  vebcat("Regions imported successfully", color = "funSuccess")
  
  return(region)
}
