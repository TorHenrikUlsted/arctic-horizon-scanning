source("./src/utils/utils.R")
source("./src/setup/setup.R")
source("./src/hypervolume/data_acquisition/acquire_data.R")
source("./src/hypervolume/data_analysis/analyze_data.R")
source("./src/hypervolume/data_processing/process_data.R")

cat(blue("Initiating hypervolume sequence \n"))

###########################################
#                                         #
# ---------------- setup ---------------- #
#                                         #
###########################################


sp_df <- setup_sp(test = T, big_test = F)


###########################################
#                                         #
# ---------- Data Acquisition ----------- #
#                                         #
###########################################

# https://nsidc.org/data/glims/data

region <- acquire_region(
  shapefiles = c(
    cavm = "./resources/region/cavm2003/cavm.shp",
    cavm_noice = "./outputs/setup/region/cavm_edited.shp"
    #glims = "./resources/region/glims_db/glims_polygons.shp"
    #rgi = "./resources/region/rgi/RGI2000-v7.0-C-05_greenland_periphery.shp",
    #rgi_sj = "./resources/region/rgi/sval_jan/RGI2000-v7.0-C-07_svalbard_jan_mayen.shp"
    #wwfEcoRegion = "./resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp"
  )
)
plot(region$cavm, col = "blue")
plot(region$cavm_noice, col = "red", add = T)
region$cavm
region$cavm <- reproject_region(region$cavm, projection = "longlat", line_issue = T)
region$glims <- reproject_region(region$glims, projection = "longlat")
identical(crs(region$cavm), crs(region$glims))

any(!is.valid(region$glims))
valid_glims <- makeValid(region$glims)
any(!is.valid(valid_glims))

glims_cavm <- crop(valid_glims, ext(region$cavm))
plot(glims_cavm)


cavm_noice <- erase(region$cavm, glims_cavm)
plot(cavm_noice)


#split test
ice2 <- terra::split(region$glims, as.factor(unlist(region$glims[[1]])))
ice3 <- terra::split(ice2[[2]], as.factor(unlist(ice2[[2]][[23]])))

glims_geo_names <- unique(ice2[[2]][[23]][[1]])
glims_geo_names <- glims_geo_names[]
for (i in seq_along(ice3)) {
  names(ice3)[i] <-  paste0(glims_geo_names[[i]])
  cat("renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(glims_geo_names[[i]]), "\n")
}

ice_cavm_list <- ice3

for (i in seq_along(ice_cavm_list)) {
  if (any(!is.valid(ice_cavm_list[[i]]))) {
    cat(red(names(ice_cavm_list[i]), "is invalid. \n"))
    validShape <- makeValid(ice_cavm_list[[i]])
    
    if (any(!is.valid(validShape))) cat(red("Failed to fix shape. \n")) else cat(green("Successfully fixed shape. \n"))
    
    ice_cavm_list[[i]] <- validShape
  } else {
    cat(green(names(ice_cavm_list[i]), "is valid. \n"))
  }
}

plot(glims_cavm)
names(ice_cavm_list)       
length(ice_cavm_list)

cavm_wo_ice <- terra::crop(ice_cavm_list$`Prince William Sound`, ext(region$cavm))
cavm_wo_ice <- crop(region$cavm, cavm_wo_ice)

# Split cavm into florisitc regions
cavm_floreg <- terra::split(region$cavm, region$cavm$FLOREG)

for (i in seq_along(cavm_floreg)) {
  names(cavm_floreg)[i] <-  paste("floreg", unique(cavm_floreg[[i]]$FLOREG), sep="_")
  cat("renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(names(cavm_floreg)[i]), "\n")
}

plot(ice_cavm_list$`Prince William Sound`)

plot(cavm_floreg$floreg_0, col = "blue")
plot(region$rgi, add = T, col = "red")


biovars_world <- acquire_biovars()

biovars_world <- scale_biovars(biovars_world)

biovars_region <- acquire_region_data(biovars_world, region, projection = "longlat")

biovars_floreg <- acquire_region_data(biovars_world, cavm_floreg, projection = "longlat")

plot(biovars_region[[1]])

acq_sp_data <- acquire_species_data(sp_df, test = T, big_test = F)

###########################################
#                                         #
# ------------ Data Analysis ------------ #
#                                         #
###########################################


analyzed_data <- analyze_correlation(biovars_region, threshold = 0.5)

# choose wanted correlation dimensions
biovars_world <- terra::subset(biovars_world, c(18, 10, 3, 4))
biovars_region <- terra::subset(biovars_region, c(18, 10, 3, 4))
for (i in seq_along(biovars_floreg)) {
  subset <- biovars_floreg[[i]][[c(18, 10, 3, 4)]]
  
  biovars_floreg[[i]] <- subset

}

###########################################
#                                         #
# ----------- Data Processing ----------- #
#                                         #
###########################################


data_processing <- function(sp_data, biovars_world, projection = "longlat", verbose = F) {
  
  cat(blue("Initiating data processing protocol \n"))
  
  tryCatch({
    cat("Processing species data \n")
    proc_sp_data <- process_sp_data(sp_data, projection = projection, verbose = verbose)
  }, error = function(e) {
    cat(red(e), "\n")
  })
  
  tryCatch({
    cat("Processing environment data \n")
    proc_env_data <- process_env_data(biovars_world, proc_sp_data, projection = projection, verbose = verbose)
  }, error = function(e) {
    cat(red(e), "\n")
  })
  
  cat("Processed environment data sample: \n")
  print(head(proc_env_data[[1]], 1))
  
  cat(cc$lightGreen("Data processing protocol completed successfully. \n"))
  
  return(proc_env_data)
}

###########################################
#                                         #
# ------------- Hypervolume ------------- #
#                                         #
###########################################

region_hv <- analyze_region_hv(biovars_region, "cavm", method = "gaussian", samples.per.point = 1, verbose = T)

hv_analysis <- function(region_hv, sp_data, method, verbose) {
  cat(blue("Initiating hypervolume sequence \n"))
  
  ## Inclusion test to eliminate obvious non-overlaps
  cat("Computing inclusion analysis. \n")
  tryCatch({
    included_sp <- analyze_inclusion(region_hv, sp_data, verbose)
    cat("Number of TRUE / FALSE values:", cc$lightSteelBlue(sum(included_sp == T)), "/", cc$lightSteelBlue(sum(included_sp == F)), "=", cc$lightSteelBlue(sum(included_sp == T) / sum(included_sp == F)))
    if(any(included_sp == T)) {
      cat(green("Included for further hypervolume analysis. \n")) 
    } else { 
      cat(red("Excluded from further hypervolume analysis. \n")) 
      return() 
    }
    
  }, error = function(e) {
    cat(red(e))
  })
  
 
  
  ## If included, continue with hypervolume analysis
  cat("Computing hypervolume analysis. \n")
  tryCatch({
    sp_hv <- hypervolume(sp_data, name = names(sp_data)[1], method = method, verbose = verbose)
    
    analyzed_hv <- analyze_hv(region_hv, sp_hv, verbose)
    
    if (analyzed_hv > 0) {
      cat(sp_hv@Name, green("Included for further hypervolume analysis. \n")) 
    } else {
      cat(sp_hv@Name, red("Excluded from further hypervolume analysis. \n")) 
      return()
    }
    
  }, error = function(e) {
    cat(red(e))
  })

  
  tryCatch({
    # Move onto projections
    create_dir_if("./outputs/hypervolume/projections")
    
    cat("Computing projection analysis. \n")
    
    ## Projections for probability
    sp_prob_project <- hypervolume_project(sp_hv, biovars_cavm, type = "probability", verbose = verbose)
    
    names(sp_prob_project) <- ("suitabilityScore")
    
    writeRaster(sp_prob_project, paste0("./outputs/hypervolume/projections/", sp_hv@Name, "/probability_longlat.tif"), format = GTiff)
    
    laea_prob_proj <- terra::project(sp_prob_project, crs(laea_crs))
    
    writeRaster(laea_prob_proj, paste0("./outputs/hypervolume/projections/", sp_hv@Name, "/probability_laea.tif"), format = GTiff)
    
    ## Projections for inclusion
    sp_inc_project <- hypervolume_project(sp_hv, biovars_cavm, type = "inclusion", verbose = verbose)
    
    names(sp_inc_project) <- ("suitabilityScore")
    
    writeRaster(sp_prob_project, paste0("./outputs/hypervolume/projections/", sp_hv@Name, "/inclusion_longlat.tif"), format = GTiff)
    
    laea_inc_proj <- terra::project(sp_inc_project, crs(laea_crs))
    
    writeRaster(laea_inc_proj, paste0("./outputs/hypervolume/projections/", sp_hv@Name, "/inclusion_laea.tif"), format = GTiff)
  }, error = function(e) {
    cat(red(e))
  })

  
  cat(cc$lightGreen("Hypervolume sequence completed successfully. \n"))
}

# joined_hv <- do.call(hypervolume_join, hv_list)


###########################################
#                                         #
# ------------- Run process ------------- #
#                                         #
###########################################

hv_log <- "./outputs/hypervolume/logs/run_log.txt"
last_iteration <- "./outputs/hypervolume/logs/last_iteration.txt"

if(file.exists(last_iteration)) {
  last_iteration <- as.integer(readLines(last_iteration))
  cat("File found, continuing from previous iteration:",last_iteration, "\n")
} else {
  last_iteration <- 0
}

for (i in (last_iteration+1):length(acq_sp_data)) {
  if (!file.exists(hv_log)) sink(hv_log) else sink(hv_log, append = T)
  
  
  
  cat("Run iteration", cc$lightSteelBlue(i), "\n")
  cat("Using species:", cc$lightSteelBlue(sp_data[[1]][[1]]), "\n")
  
  
  processed_data <- data_processing(acq_sp_data[[i]], biovars_world)
  analyzed_hv <- hv_analysis(region_hv, processed_data, method, verbose)
  
  # Append to csv file
  create_dir_if("")
  fwrite("./outputs/hypervolume/", append = T, bom = T)
  
  sink()
  
  writeLines(as.character(i), last_iteration)
}
