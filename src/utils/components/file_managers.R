create_dir_if <- function(dirs, keep = TRUE) {
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
    } else {
      if (keep == FALSE) {
        unlink(d, recursive = TRUE)

        dir.create(d, recursive = TRUE)
      }
    }
  }
}

create_file_if <- function(files, keep = F) {
  for (f in files) {
    if (!file.exists(f)) {
      file.create(f)
    } else {
      if (keep == FALSE) {
        file.remove(f)
        
        file.create(f)
      }
    }
  }
}

download_region_if <- function(region.file, download.link, download.page) {
  
  file_ext <- tail(strsplit(region.file, split = "\\.")[[1]], 1)
  create_dir_if(dirname(region.file))
  
  if (!file.exists(region.file)) {
    vebcat("Could not find region.", color = "nonFatalError")
    catn("Downloading...")
    
    download_status <- try(download.file(download.link, destfile = region.file), silent = TRUE)
    
    if (inherits(download_status, "try-error")) {
      vebcat("Opening download page for manual download.", color = "indicator")
      
      Sys.sleep(1000)
      
      catn("Upload the file to:", colcat(region.file, color = "indicator"))
      
      Sys.sleep(2000)
      
      utils::browseURL(download.page)
    }
    
    if (file_ext == "zip") {
      catn("File downloaded, unzipping...")
      unzip(region.file, exdir = dirname(region.file))
    }
  }
}
