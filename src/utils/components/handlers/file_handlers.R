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
  vebcat("Detecting installation of", name, color = "funInit")
  
  file_ext <- tail(strsplit(basename(out.file), split = "\\.")[[1]], 1)
  create_dir_if(dirname(out.file))
  
  vebprint(name, verbose, "Name of file:")
  vebprint(file_ext, verbose, "Extension of file:")
  
  file_save <- paste0(dirname(out.file), "/", name, ".", download.file.ext)
  save_ext <- tail(strsplit(basename(file_save), split = "\\.")[[1]], 1)
  
  vebprint(file_save, verbose, "Name of file to be downloaded:")
  vebprint(save_ext, verbose, "Extension of file to be downloaded:")
  
  if (!file.exists(out.file)) {
    if(!is.null(download.direct)) {
      catn("Installation not found, Downloading...")
      
      download_status <- try(download.file(download.direct, destfile = file_save, mode = "wb", quiet = TRUE))
      
      if (save_ext == "zip") {
        catn("File downloaded, unzipping...")
        unzip(file_save, exdir = dirname(out.file))
        
        # Only remove zip file AFTER confirming the extracted files exist
        if (file.exists(out.file)) {
          file.remove(file_save)
        } else {
          stop("Failed to extract files from zip archive")
        }
        
        # Handle nested directory structure
        new_dir <- list.files(dirname(out.file), full.names = TRUE, include.dirs = TRUE)
        new_dir <- new_dir[file.info(new_dir)$isdir]  # Get only directories
        
        # Check for new directory creation
        if (length(new_dir) > 0) {
          vebprint(new_dir, verbose, text = "Renaming files from new dir:")
          vebprint(dirname(out.file), verbose, text = "To out dir:")
          catn("Files unzipped, renaming...")
          
          # Move all files from subdirectories
          for (dir in new_dir) {
            files <- list.files(dir, full.names = TRUE, recursive = TRUE)
            for (file in files) {
              filename <- basename(file)
              out_file <- file.path(dirname(out.file), filename)
              file.copy(from = file, to = out_file, overwrite = TRUE)
            }
          }
          
          catn("Removing extra directories.")
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
    vebcat(name, "Already installed at", out.file, color = "funSuccess")
  }
}

import_font_if <- function(font, paths, pattern) {
  tryCatch({
    # Check if font already exists
    if (font %in% fonts()) {
      catn(font, "already imported")
      return(invisible())
    }
    
    # Create directory if it doesn't exist
    if (!dir.exists(paths)) {
      dir.create(paths, recursive = TRUE)
    }
    
    # List font files
    font_files <- list.files(paths, pattern = pattern, full.names = TRUE)
    
    # If no font files found, report error
    if (length(font_files) == 0) {
      stop("No font files found in ", paths, " matching pattern ", pattern)
    }
    
    # Import fonts
    catn("Importing font:", font)
    font_import(
      paths = paths,
      pattern = pattern,
      prompt = FALSE  # Don't prompt for confirmation
    )
    
    # Reload font database
    loadfonts(device = "all", quiet = TRUE)
    
  }, error = function(e) {
    warning("Font import failed: ", e$message, "\n",
            "Continuing without custom font.")
  })
}

download_github_dir_if <- function(repo.owner, repo.name, branch = "main", dir.path, dir.out, file.exclude = NULL, verbose = FALSE) {
  
  vebcat("Detecting", repo.name, "installation...", color = "funInit")
  
  dir_out_len <- length(list.files(dir.out))
  
  if (dir.exists(dir.out) &  dir_out_len > 0) {
    catn("Location", colcat(dir.out, color = "output"))
    vebcat("GitHub repository", highcat(repo.name), "already downloaded with", highcat(dir_out_len), "files.", color = "funSuccess")
    return(invisible())
  }
  
  catn("Installation not found, downloading...")
  
  create_dir_if(dir.out)
  
  # Construct the API URL
  api_url <- paste0("https://api.github.com/repos/", repo.owner, "/", repo.name, "/contents/", dir.path, "?ref=", branch)
  
  vebprint(api_url, verbose, "API URL:")
  
  # Make the API request
  response <- GET(api_url)
  
  # Check if the request was successful
  if (status_code(response) != 200) {
    stop("Failed to retrieve folder contents. Status code: ", status_code(response))
  }
  
  # Parse the JSON response
  contents <- data.table(fromJSON(rawToChar(response$content)))
  
  exclude_pattern <- paste(file.exclude, collapse = "|")
  
  if (!is.null(file.exclude)) {
    contents <- contents[!grepl(exclude_pattern, contents$name, ignore.case = TRUE)]
  }
  
  vebprint(contents, verbose, "Contents:")
  
  # Download each file
  for (i in 1:nrow(contents)) {
    item <- contents[i]
    if (item$type == "file") {
      item_name <- item$name
      catn("Downloading", item_name)
      download_url <- item$download_url
      out_path <- paste0(dir.out, "/", item$name)
      download.file(download_url, out_path, mode = "wb")
    }
  }
  
  catn("All GitHub files have been downloaded to:", colcat(dir.out, color = "output"))
  vebcat("Successfully installed", repo.name, color = "funSuccess")
}
