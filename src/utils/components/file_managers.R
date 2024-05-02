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
      if (!dir.exists(dirname(f))) {
        dir.create(dirname(f), recursive = TRUE)
      }
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
  region_name <- basename(region.file)
  vebcat("Checking need for download for", region_name, color = "funInit")
  
  file_ext <- tail(strsplit(region.file, split = "\\.")[[1]], 1)
  create_dir_if(dirname(region.file))
  
  if (!file.exists(region.file)) {
    vebcat("Could not find region.", color = "nonFatalError")
    
    if(!is.null(download.link)) {
      catn("Downloading...")
      
      download_status <- try(download.file(download.link, destfile = region.file), silent = TRUE)
      
      if (file_ext == "zip") {
        catn("File downloaded, unzipping...")
        unzip(region.file, exdir = dirname(region.file))
      }
      
      if (inherits(download_status, "try-error")) {
        vebcat("Opening download page for manual download.", color = "indicator")
        
        Sys.sleep(1)
        
        catn("Upload the file(s) to:", colcat(dirname(region.file), color = "indicator"))
        
        Sys.sleep(2)
        
        utils::browseURL(download.page)
        
        stop()
      }
    } else {
      vebcat("Opening download page for manual download.", color = "indicator")
      
      Sys.sleep(1)
      
      catn("Upload the file(s) to:", colcat(dirname(region.file), color = "indicator"))
      
      Sys.sleep(2)
      
      utils::browseURL(download.page)
      
      stop("Could not automatically download region. Stopping...")
    }
    
  } else {
    vebcat("No need to download", region_name, color = "funSuccess")
  }
}
