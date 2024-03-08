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
