if(!file.exists("resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp")) {
  downld_url = "https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip"
  fallbck_url = "https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world"
  paper_url = "https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2"
  zip_file = "resources/region/wwfTerrestrialEcoRegions.zip"
  dest_dir = "resources/region/wwfTerrestrialEcoRegions"
  official_dir = file.path("resources/region/wwfTerrestrialEcoRegions", "official")
  
  
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