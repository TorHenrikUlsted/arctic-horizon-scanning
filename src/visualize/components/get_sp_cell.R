get_sp_cell <- function(rast, region, method, cores.max = 1, prob.threshold = NULL) {
  cat(blue("Acquiring number of species per cell.\n"))
  
  # Define the function to apply to each chunk
  if (verbose) cat("Using app() to calculate included and excluded species in each cell. \n")
  
  if (method == "inclusion") {
    sp_count <- terra::app(rast, fun=function(x) {
      included = sum(!is.na(x) & x == 1)
      excluded = sum(!is.na(x) & x == 0)
      prop = ifelse(included == 0 & excluded == 0, 0, included / (included + excluded))
      
      c(included = included, excluded = excluded, prop = prop)
    })
  } else if (method == "probability") {
      sp_count <- terra::app(rast, fun=function(x) {
        if (!is.null(prob.threshold)){
          prob_max = max(!is.na(x) & x > prob.threshold)
          prob_mean = mean(!is.na(x) & x > prob.threshold)
          prob_median = median(!is.na(x) & x > prob.threshold)
        } else {
          prob_max = max(!is.na(x) & x)
          prob_mean = mean(!is.na(x) & x)
          prob_median = median(!is.na(x) & x)
        } 
          c(prob_max = prob_max, prob_mean = prob_mean, prob_median = prob_median)
      })
  }
  
  
  sp_count_dt <- data.table(terra::values(sp_count, na.rm = TRUE))
  
  # Remove rows where both included and excluded are 0
  if (verbose) cat("Removing rows where all are 0. \n")
  if (method == "inclusion") {
    sp_count_dt <- sp_count_dt[!(included == 0 & excluded == 0), ]
  } else if (method == "probability") {
    sp_count_dt <- sp_count_dt[!(prob_max == 0 & prob_mean == 0 & prob_median == 0), ]
  }
  
  if (verbose) cat("Append lon/lat values to region.\n")
  coords_dt <- data.table(terra::xyFromCell(rast, 1:nrow(sp_count_dt)))
  setnames(coords_dt, c("lon", "lat"))
  sp_count_dt <- cbind(sp_count_dt, coords_dt)

  if (method == "inclusion") {
    if (verbose) cat("Calculating proportion of included species. \n")
    sp_count_dt$incProp <- sp_count_dt$included / (sp_count_dt$included + sp_count_dt$excluded)
  }

  ####################
  #     Region       #
  ####################

  if (verbose) cat("Extracting cells for the different regions. \n")
  region_cell_vals <- terra::extract(sp_count, region)
  
  region_cell_vals <- as.data.table(region_cell_vals, na.rm = TRUE)
  if (method == "inclusion") {
    colnames(region_cell_vals) <- c("regionId", "included", "excluded", "prop")
  } else if (method == "probability") {
    print(region_cell_vals)
    colnames(region_cell_vals) <- c("regionId", "prob_max", "prob_mean", "prob_median")
  }
  
  region_dt <- as.data.table(region)
  
  # Add a region id to the cavm
  if (verbose) cat("Adding regionId to the SpatVector. \n")
  region_dt[, regionId := .I]
  
  if (verbose) cat("Merging data tables by regionId. \n")
  sub_regions_dt <- merge(region_cell_vals, region_dt, by = "regionId")
  
  if (verbose) cat("Removing rows where all are 0. \n")

  if (method == "inclusion") {
    sub_regions_dt <- sub_regions_dt[!(included == 0 & excluded == 0), ]
  } else if (method == "probability") {
    sub_regions_dt <- sub_regions_dt[!(prob_max == 0 & prob_mean == 0 & prob_median == 0), ]
  }
  
  if (verbose) cat("Append lon/lat values to sub-region.\n")
  coords_dt <- data.table(terra::xyFromCell(rast, 1:nrow(sub_regions_dt)))
  setnames(coords_dt, c("lon", "lat"))
  sub_regions_dt <- cbind(sub_regions_dt, coords_dt)
  
  if (method == "inclusion") {
    if (verbose) cat("Calculating proportion of included species. \n")
    sub_regions_dt$incProp <- sub_regions_dt$included / (sub_regions_dt$included + sub_regions_dt$excluded)
  }
  
  cat(cc$lightGreen("Successfully calculated species in each cell. \n"))
  
  return(list(
    sp_cell_arctic = sp_count_dt,
    sp_count = sp_count,
    sp_cell_floreg = sub_regions_dt,
    region_count = region_cell_vals
  ))
}