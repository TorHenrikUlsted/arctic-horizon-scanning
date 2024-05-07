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

download_if <- function(out.file, download.file.ext, download.direct, download.page) {
  name <- sub("\\..*$", "", basename(out.file))
  vebcat("Checking need for download for", name, color = "funInit")
  
  file_ext <- tail(strsplit(basename(out.file), split = "\\.")[[1]], 1)
  create_dir_if(dirname(out.file))
  
  file_save <- gsub(file_ext, download.file.ext, out.file)
  save_ext <- tail(strsplit(basename(file_save), split = "\\.")[[1]], 1)
  
  if (!file.exists(out.file)) {
    if(!is.null(download.direct)) {
      catn("Downloading...")
      
      download_status <- try(download.file(download.direct, destfile = file_save, silent = TRUE))
      
      if (save_ext == "zip") {
        catn("File downloaded, unzipping...")
        unzip(file_save, exdir = dirname(out.file))
        file.remove(file_save)
        
        new_dir <- paste0(dirname(out.file), "/", list.files(dirname(out.file)))
        
        for (file in list.files(new_dir, full.names = TRUE)) {
          filename <- basename(file)
          out_file <- paste0(dirname(out.file), "/", filename)
          file.rename(from = file, to = out_file)
        }
        
        unlink(new_dir, recursive = TRUE)
      }
      
      if (inherits(download_status, "try-error")) {
        vebcat("Opening download page for manual download.", color = "indicator")
        
        Sys.sleep(1)
        
        catn("Upload the file(s) to:", colcat(dirname(out.file), color = "indicator"))
        
        Sys.sleep(2)
        
        utils::browseURL(download.page)
        
        stop()
      }
    } else {
      vebcat("Opening download page for manual download.", color = "indicator")
      
      Sys.sleep(1)
      
      catn("Upload the file(s) to:", colcat(dirname(out.file), color = "indicator"))
      
      Sys.sleep(2)
      
      utils::browseURL(download.page)
      
      stop("Could not automatically download file. Stopping...")
    }
    
  } else {
    vebcat("No need to download", name, color = "funSuccess")
  }
}
