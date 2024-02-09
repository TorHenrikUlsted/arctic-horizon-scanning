check_hv_output <- function(spec.list, hv.dir, hv.method, hv.projection = laea_crs, hv.inc.t = 0.5, hv.clean = F, verbose = F) {
  
  proj_dir <- paste0(hv.dir, "/projections/", hv.method)
  stats_dir <- paste0(hv.dir, "/stats")
  logs_dir <- "./outputs/visualize/logs"
  missing_rast_log <- paste0(logs_dir, "/missing-rasters.txt")
  missing_sp_log <- paste0(logs_dir, "/missing-species.txt")
  
  sp_dirs <- list.dirs(path = proj_dir, full.names = TRUE)
  sp_dirs_dt <- data.table(species = basename(sp_dirs))
  
  spec_list_dt <- data.table(dir_name = basename(spec.list))
  spec_list_dt$iteration <- seq_along(spec_list_dt)
  sp_it_num <- c()
  
  
  for (i in seq_along(sp_dirs)) {
    sp_dir <- sp_dirs[[i]]
    
    rast_files <- list.files(path = sp_dir, pattern = "\\.tif$", full.names = TRUE)
    
    
    # If the dir has less files than expected, add the iteration of the species to a file
    if (length(rast_files) < length(hv.inc.t) + 1) {
      sp_it_num <- spec_list_dt$iteration[match(basename(sp_dir), spec_list_dt$dir_name)]
      
      try(mr_con <- file(missing_rast_log, open = "at"))
        writeLines(as.character(sp_it_num), con = mr_con)
      close(mr_con)
    }
    
    cat("\rChecking directory", i, "/", length(sp_dirs))
    flush.console()
    
    for (j in seq_along(rast_files)) {
      rast_file <- rast_files[[j]]
      
      raster_crs <- terra::crs(rast_file)
      
      if (raster_crs != hv.projection && hv.clean) {
        if (verbose) cat("Raster is in the wrong projection. \n")
        raster <- terra::rast(rast_file)
        raster <- terra::project(raster, projection.expected)
        terra::writeRaster(raster, rast_file)
      } else {
        if (verbose) cat("Raster is in the correct projection. \n")
      }
    }
    
  }
  
  # Check if stats has the same length as the number of species directories
  stat_files <- list.files(stats_dir, full.names = TRUE)
    
  for (i in seq_along(stat_files)) {
    stat <- stat_files[[i]]
    
    if (verbose) {cat("stat:\n"); print(stat)}
    
    dt <- fread(stat)
    
    if (nrow(dt) != length(sp_dirs)) {
      cat(cc$lightCoral(basename(stat), "failed the length check.\n"))
      cat("    Length of stats / length of species:\n")
      cat("   ",cc$lightSteelBlue(nrow(dt)), "/", cc$lightSteelBlue(length(sp_dirs)), "\n")
      
      missing_sp <- setdiff(stat$species, sp_dirs_dt$species)
      
      if (verbose) {cat("matched_sp:\n"); print(matched_sp)}
      
      if (length(missing_files) > 0) {
        # Get the iteration number from spec_list
        missing_files_dt <- data.table(dir_name = missing_files)
        missing_files_dt <- merge(missing_files_dt, spec_list_dt, by = "dir_name")
        
        try(msp_con <- file(missing_sp_log, open = "at"))
          writeLines(as.character(missing_files_dt$iteration), con = msp_con)
        close(msp_con)
        
      } else {
        cat(cc$lightGreen(basename(stat), "passed the length check.\n"))
      }
      
    }
    
  }

}