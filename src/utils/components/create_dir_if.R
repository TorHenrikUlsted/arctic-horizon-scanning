create_dir_if <- function(d, keep = T) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = T)
  } else {
    if (keep == F) {
      unlink(d, recursive = T)
      dir.create(d, recursive = T)
    }
  }
}