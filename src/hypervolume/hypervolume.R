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
    glims = "./resources/region/glims_db/glims_polygons.shp",
    rgi = "./resources/region/rgi/RGI2000-v7.0-C-05_greenland_periphery.shp",
    rgi_sj = "./resources/region/rgi/sval_jan/RGI2000-v7.0-C-07_svalbard_jan_mayen.shp"
    #wwfEcoRegion = "./resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp"
  )
)

region$cavm <- reproject_region(region$cavm, projection = "laea", line_issue = T)
plot(region$rgi_sj)
plot(region$cavm)
plot(region$rgi, add = T, col = "red")

cavm_wo_ice_green <- erase(region$cavm, region$rgi)
plot(cavm_wo_ice_green, col = "red")

ice <- terra::project(region$glims, crs(region$cavm))
identical(crs(region$cavm), crs(ice))
ice2 <- terra::split(region$glims, as.factor(unlist(region$glims[[1]])))

names(ice2[[2]])
unique(ice2[[2]][[23]][[1]][[1]])

ice3 <- terra::split(ice2[[2]], as.factor(unlist(ice2[[2]][[23]])))
plot(ice3[[33]])

glims_geo_names <- unique(ice2[[2]][[23]][[1]])
glims_geo_names <- glims_geo_names[]
for (i in seq_along(ice3)) {
  names(ice3)[i] <-  paste0(glims_geo_names[[i]])
  cat("renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(glims_geo_names[[i]]), "\n")
}

plot(ice3$`Prince William Sound`)

terra::intersect(region$cavm, ice3[[2]])
ice_cavm_list <- ice3

for (i in seq_along(ice3)) {
    
  
  if(length(intersect(ice3[[i]], region$cavm)) == 0) {
    
    ice_cavm_list[[i]] <- NULL
  }

}

# Remove ice covered areas from the CAVM
cavm_wo_ice <- terra::crop(region$cavm, region$glims)

cavm_floreg <- terra::split(region$cavm, region$cavm$FLOREG)

for (i in seq_along(cavm_floreg)) {
  names(cavm_floreg)[i] <-  paste("floreg", unique(cavm_floreg[[i]]$FLOREG), sep="_")
  cat("renaming item", cc$lightSteelBlue(i), "to", cc$lightSteelBlue(names(cavm_floreg)[i]), "\n")
}

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

proc_sp_data <- process_sp_data(acq_sp_data, projection = "longlat", verbose = T)

proc_env_data <- process_env_data(biovars_world, proc_sp_data, projection = "longlat", verbose = T)

cat("Processed environment data sample: \n")
print(head(proc_env_data[[1]], 3))


###########################################
#                                         #
# ------------- Hypervolume ------------- #
#                                         #
###########################################

cat("Running hypervolume for", cc$lightSteelBlue(names(proc_env_data)[1])) 
cat("Samples per point", cc$lightSteelBlue(ceiling((10^(3 + sqrt(ncol(proc_env_data[[1]]))))/nrow(proc_env_data[[1]]))), "\n")


## Inclusion test to eliminate obvious non-overlaps
inc_test <- hypervolume_inclusion_test(region_hv_box, proc_env_data[[1]], reduction.factor = 1, fast.or.accurate =
                                         "accurate", fast.method.distance.factor = 1,
                                       accurate.method.threshold =
                                         quantile(region_hv_box@ValueAtRandomPoints,
                                                  0.5), verbose = TRUE)

cat("Number of TRUE / FALSE values:", cc$lightSteelBlue(sum(inc_test == T)), "/", cc$lightSteelBlue(sum(inc_test == F)), "=", cc$lightSteelBlue(sum(inc_test == T) / sum(inc_test == F)))
if(any(inc_test == T)) cat(green("Species added to further hypervolume testing. \n")) else cat(red("Species is removed from further testing. \n"))



## box probability

sp_hv_box <- hypervolume_box(proc_env_data[[1]], name = "species", verbose = T)
region_hv_box <- analyze_region_hv(biovars_region, "cavm", method = "box", samples.per.point = 1, verbose = T)
hv_set <- hypervolume_set(sp_hv_box, region_hv, check.memory=F)
hv_stats <- hypervolume_overlap_statistics(hv_set)

sp_surviv_region <- 1 - hv_stats[[4]]


svalbard_hv <- hypervolume_box(biovars_floreg[[11]], name = "svalbard", verbose = T)
sval_set <- hypervolume_set(sp_hv_box, svalbard_hv, check.memory=F)
hypervolume_overlap_statistics(sval_set)



draw.pairwise.venn(
  area1 = hv_stats[[3]] + hv_stats[[1]],
  area2 = hv_stats[[4]] + hv_stats[[1]],
  cross.area =  hv_stats[[1]],
  category = c(sp_hv_box@Name, region_hv_box@Name),
  fill = c("blue", "red"),
  lty = "blank"
)

floreg_project <- hypervolume_project(sp_hv_box, biovars_floreg[[11]], type = "probability", verbose = T)
plot(floreg_project)

plot(biovars_floreg[[11]][[1]])
sp_hv_gauss <- hypervolume_gaussian(proc_env_data[[1]], name = "species", quantile.requested.type = "probability", verbose = T)

floreg_gauss_proj <- hypervolume_project(sp_hv_gauss, biovars_floreg[[1]], type = "probability", verbose = T)
plot(floreg_gauss_proj)

fgps <- hypervolume_project(sp_hv_gauss, biovars_floreg[[11]], type = "probability", verbose = T)
plot(fgps)

hv_set <- hypervolume_set(sp_hv_gauss, region_hv, check.memory=F)
hv_stats <- hypervolume_overlap_statistics(hv_set)

sp_hvs <- analyze_hypervolume(proc_env_data, sample_size = 1, verbose = T)

plot_hypervolumes(sp_hvs@HVList)
plot(sp_hvs@HVList[[1]])

region_hv <- analyze_region_hv(biovars_region, "cavm", method = "box", samples.per.point = 1, verbose = T)


png("projection_gaussian_svaldbard.png", width = 1000, height = 500)
plot(fgps, legend = "right", box = F)
title(main = "Probability density of Saxifraga in CAVM", outer = F)
dev.off()

bt <- terra::project(hv_project, crs("+proj=laea +lon_0=180 +lat_0=90 +datum=WGS84"))
plot(bt, legend = "right", box = F)

