get_sp_cell <- function(rast, shape, method, cores.max = 1, prob.threshold = NULL, out.filename, verbose = F) {
  cat(blue("Acquiring number of species per cell.\n"))
  
  # Init files and folders
  create_dir_if(dirname(out.filename))
  
  if (method == "inclusion") {
    if (!file.exists(out.filename)) {
      if (verbose) {
        cat("out.filename:", out.filename, "\n")
        cat("rast:\n")
        print(rast)
      }

      cat("Using app() to calculate included and excluded species in each cell. \n")      
      sp_count <- terra::app(rast, cores = cores.max, fun=function(x) { ## Cores Ignores C++ function commands like sum(), mean() etc..
        included = sum(!is.na(x) & x > 0)
        excluded = length(x) - ceiling(included)
        prop = ifelse(included == 0 & excluded == 0, 0, included / (included + excluded))
        
        c(included = included, excluded = excluded, prop = prop)
      })
      
      cat("Writing out to raster.\n")
      writeRaster(sp_count, out.filename, overwrite = TRUE)
    } else {
      cat(out.filename, "found, reading file.\n")
      sp_count <- rast(out.filename)
    }

  } else if (method == "probability") {
    if(!file.exists(out.filename)) {
      cat("Using app() to calculate probability for species in each cell. \n")
      sp_count <- terra::app(rast, cores = cores.max, fun=function(x) {
        if (!is.null(prob.threshold)){
          total = sum(x[!is.na(x) & x > prob.threshold])
          prop = ifelse(total == 0, 0, x / total)
        } else {
          total = sum(x[!is.na(x)])
          prop = ifelse(total == 0, 0, x / total)
        } 
        c(total = total, prop = prop)
      })
      
      cat("Writing to raster. \n")
      
      writeRaster(sp_count, out.filename, overwrite = TRUE)
    }
    cat(out.filename, "found, reading file.\n")
    sp_count <- rast(out.filename)
  }
  
  region_dt <- data.table(terra::values(sp_count, na.rm = TRUE))

  ####################
  #     Region       #
  ####################


  
  
  return(list(
    region = sub_regions_dt,
    sp_region = sp_regions_dt,
    sp_count = sp_count
  ))
}

extract_region <- function(sp_count, shape) {
  cat("Extracting cells for the different regions. \n")
  region_cell_vals <- terra::extract(sp_count, shape, cells = TRUE)
  
  region_cell_vals <- as.data.table(region_cell_vals, na.rm = TRUE)
  
  if (method == "inclusion") {
    colnames(region_cell_vals) <- c("regionId", "included", "excluded", "prop", "cell")
    
  } else if (method == "probability") {
    print(region_cell_vals)
    colnames(region_cell_vals) <- c("regionId", "total", "prop", "cell")
  }
  
  if (verbose)  {
    cat("region_cell_vals:\n")
    print(region_cell_vals)
  }
  
  shape_dt <- as.data.table(shape)
  
  # Add a region id to the cavm
  if (verbose) cat("Adding regionId to the SpatVector. \n")
  shape_dt[, regionId := .I]
  
  if (verbose) {
    cat("shape_dt:\n")
    print(shape_dt)
  }
  
  if (verbose) cat("Merging data tables by regionId. \n")
  sub_regions_dt <- merge(region_cell_vals, shape_dt[, .(regionId, FLOREG, country, floristicProvince)], by = "regionId")
  
  if (verbose) {
    cat("sub_regions_dt:\n")
    print(sub_regions_dt)
  }
  
}

extract_species <- function() {
  cat("Extracting species names from each cell. \n")
  sp_name_cell <- terra::extract(rast, shape, cells=TRUE)
  
  print(sp_name_cell)
  
  sp_dt <- as.data.table(sp_name_cell)
  
  # #Since ID is a column it will be in the long format, remove it
  sp_cols <- setdiff(names(sp_dt), c("ID"))  # Exclude 'ID' and 'cell'
  sp_dt <- sp_dt[, ..sp_cols]
  sp_dt[, cell := .I]
  
  # Append the cell ID and reshape
  
  if (verbose)  {
    cat("sp_dt:\n")
    print(sp_dt)
  }
  
  sp_dt_long <- melt(sp_dt, id.vars = "cell", variable.name = "species", value.name = "presence")
  sp_dt_long <- sp_dt_long[presence > 0, .(cell, species)]
  
  if (verbose) {
    cat("melted species:\n")
    print(sp_dt_long)
  }
  
  sp_regions_dt <- merge(sub_regions_dt, sp_dt_long, by = "cell")
  
  if (verbose) {
    cat("sp_regions_dt:\n")
    print(sp_regions_dt)
    
    cat("Append lon/lat values to sub-region.\n")
  }
  
  cat(cc$lightGreen("Successfully calculated species in each cell. \n"))
}