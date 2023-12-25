get_disk_space <- function(directory = "/", units = "B") {
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    disk_space <- system2("cmd", args = "/c fsutil volume diskfree C:", stdout = TRUE)
    disk_space <- as.numeric(strsplit(disk_space[3], ": ")[[1]][2])
  } else if (os %in% c("Linux", "Darwin")) {
    disk_space <- system(paste0("df --output=avail -h ", directory, " | tail -n 1"), intern = TRUE)
    disk_space <- trimws(disk_space)  # Remove whitespace
    unit <- substr(disk_space, nchar(disk_space), nchar(disk_space))
    disk_space <- as.numeric(gsub(unit, "", disk_space))
    if (unit == "T") {
      disk_space <- disk_space * 1024  # Convert from terabytes to gigabytes
    }
  } else {
    stop("Unsupported operating system.")
  }
  
  conversion_factors <- c(B = 1, KB = 1024, MB = 1024^2, GB = 1024^3, TB = 1024^4)
  disk_space <- disk_space * 1e9 / conversion_factors[units]
  
  return(disk_space)
}
