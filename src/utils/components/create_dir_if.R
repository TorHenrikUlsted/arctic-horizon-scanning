create_dir_if <- function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = T)
}