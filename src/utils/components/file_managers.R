create_dir_if <- function(..., keep = TRUE) {
  dirs <- list(...)
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

create_file_if <- function(..., keep = FALSE) {
  files <- list(...)
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

download_if <- function(out.file, download.file.ext, download.direct = NULL, download.page, verbose = FALSE) {
  name <- sub("\\..*$", "", basename(out.file))
  vebcat("Checking need for download for", name, color = "funInit")
  
  file_ext <- tail(strsplit(basename(out.file), split = "\\.")[[1]], 1)
  create_dir_if(dirname(out.file))
  
  vebprint(name, verbose, "Name of file:")
  vebprint(file_ext, verbose, "Extension of file:")
  
  file_save <- gsub(file_ext, download.file.ext, out.file)
  save_ext <- tail(strsplit(basename(file_save), split = "\\.")[[1]], 1)
  
  vebprint(file_save, verbose, "Name of file to be downloaded:")
  vebprint(save_ext, verbose, "Extension of file to be downloaded:")
  
  if (!file.exists(out.file)) {
    if(!is.null(download.direct)) {
      catn("Downloading...")
      
      download_status <- try(download.file(download.direct, destfile = file_save, silent = TRUE))
      
      if (save_ext == "zip") {
        catn("File downloaded, unzipping...")
        unzip(file_save, exdir = dirname(out.file))
        file.remove(file_save)
        
        new_dir <- paste0(dirname(out.file), "/", list.files(dirname(out.file)))
        
        # check for new dir creation
        if (!identical(dirname(out.file), dirname(new_dir[1]))) {
          vebprint(dirname(new_dir[1]), verbose, "Renaming files to dir:")
          catn("Files unzipped, renaming...")
          
          for (file in list.files(new_dir, full.names = TRUE)) {
            filename <- basename(file)
            out_file <- paste0(dirname(out.file), "/", filename)
            file.rename(from = file, to = out_file)
          }
          
          catn("Removing extra directory.")
          
          unlink(new_dir, recursive = TRUE)
        }
        
        vebcat("File unzipped successfully", color = "proSuccess")
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

import_font_if <- function(font, paths, pattern) {
  if (font %in% fonts()) {
    catn(font, "already imported")
  } else {
    font_import(
      paths = paths,
      pattern = pattern
    )
  }
}
