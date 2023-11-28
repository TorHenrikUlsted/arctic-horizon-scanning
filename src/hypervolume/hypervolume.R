source("./src/utils/utils.R")
source("./src/setup/setup.R")
source("./src/hypervolume/data_acquisition/acquire_data.R")
source("./src/hypervolume/data_analysis/analyze_data.R")
source("./src/hypervolume/data_processing/process_data.R")

cat(blue("Initiating hypervolume sequence \n"))

sp_df <- setup_env(test = T, big_test = F)

shapefiles = c(
  cavm = "./resources/region/cavm2003/cavm.shp",
  wwfEcoRegion = "./resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp"
)

acquired_data <- acquire_data(shapefiles, projection = "laea", sp_df, plot_it = F, verbose = T)

# Rename acquired data
biovars_cavm <- acquired_data[[1]]
biovars_world <- acquired_data[[2]]
sp_data <- acquired_data[[3]]

analyzed_data <- analyze_correlation(biovars_cavm, threshold = 0.5)

# choose wanted correlation values
biovars_world <- terra::subset(biovars_world, c(18, 10, 3, 4))
biovars_cavm <- terra::subset(biovars_cavm, c(18, 10, 3, 4))

env_data <- process_data(biovars_world, sp_data, projection = "longlat", verbose = T)
  
print(head(env_data$`Saxifraga oppositifolia`, 3))

#Convert points to laea
sp_points_laea <- terra::project(env_data[[2]][[1]], crs(biovars_cavm))

plot(biovars_cavm[[2]])
plot(sp_points_laea, add =T, col = "blue", cex = .3)

plot(biovars_world[[2]])
plot(env_data[[2]][[1]], add =T, col = "blue", cex = .2)


###########################################
#                                         #
#                                         #
# ---------------- DynRB ---------------- #
#                                         #
#                                         #
###########################################


saxi_df <- as.data.frame(env_data$`Saxifraga oppositifolia`)
saxi_df$name <- "Saxifraga oppositifolia"
saxi_df <- saxi_df[, c("name", "bio_18", "bio_10", "bio_3", "bio_4")]

cavm_df <- as.data.frame(biovars_cavm)
cavm_df$name <- "cavm"
cavm_df <- cavm_df[, c("name", "bio_18", "bio_10", "bio_3", "bio_4")]

saxi_dynrb <- dynRB_VPa(rbind(saxi_df, cavm_df), steps = 201)

saxi_dynrb_df <- as.data.frame(saxi_dynrb$plot_data_niche_size_V1)

ggplot(df, aes(x = V1, y = V2)) +
  geom_point() +
  theme_minimal()

View(saxi_dynrb)

str(saxi_dynrb$plot_data_niche_size_V1)


theme_change <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  panel.border = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_text(colour="black", size = rel(1.5), angle=35, hjust = 1),
  axis.text.y = element_text(colour="black", size = rel(1.5)),
  axis.title.x = element_blank(),
  axis.title.y = element_blank()
)
result <- saxi_dynrb$result
Overlap <- as.numeric(ifelse(result$V1 == result$V2, "NA", result$port_prod))  
# 'result$port_prod' may be changed to 'result$port_mean' or 'result$port_gmean'
is.numeric(Overlap)
Result2<-cbind(result, Overlap)
breaks <- seq(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE),  
              by=round(max(Overlap, na.rm=TRUE)/10, digits=3))
col1 <- colorRampPalette(c("white", "navyblue")) #define color gradient
ggplot(Result2, aes(x = V1, y = V2)) +
  geom_tile(data = subset(Result2, !is.na(Overlap)), aes(fill = Overlap), color="black") +
  geom_tile(data = subset(Result2,  is.na(Overlap)), fill = "lightgrey", color="black") +
  scale_fill_gradientn(colours=col1(8), breaks=breaks, guide="colorbar",  
                       limits=c(min(Overlap, na.rm=TRUE),max(Overlap, na.rm=TRUE))) +
  theme_change

venn <- draw.pairwise.venn(
  area1 = result$vol_V1_gmean[[2]],
  area2 = result$vol_V2_gmean[[2]],
  cross.area = result$port_gmean[[2]],
  category = c(unique(saxi_df$name), unique(cavm_df$name)),
  fill = c("blue", "red"),
  lty = "blank"
)

grid.draw(venn)


###########################################
#                                         #
#                                         #
# ------------- Hypervolume ------------- #
#                                         #
#                                         #
###########################################


str(env_data)
print(head(env_data[[1]], 5))

ceiling((10^(3 + sqrt(ncol(env_data[[1]]))))/nrow(env_data[[1]]))
ceiling((10^(3 + sqrt(ncol(biovars_cavm))))/nrow(biovars_cavm))

estimate_bandwidth(env_data[[1]])
estimate_bandwidth(biovars_cavm)

ind_hvs <- get_ind_hv(env_data)

sp_hvs <- analyze_hypervolume(env_data, sample_size = 1000, verbose = T)

plot_hypervolumes(sp_hvs@HVList)
plot(sp_hvs@HVList[[2]])

cavm_hv <- analyze_region_hv(biovars_cavm, "cavm", verbose = T)

cavm_ind_hv <- get_ind_rs(biovars_cavm, "cavm")

plot(cavm_ind_hv@HVList[[4]])

length(ind_hvs$`Saxifraga oppositifolia`@HVList[[1]]@Data)
any(is.na(ind_hvs$`Saxifraga oppositifolia`@HVList[[1]]@Data))
length(ind_hvs$`Saxifraga oppositifolia`@HVList[[2]]@Data)
any(is.na(ind_hvs$`Saxifraga oppositifolia`@HVList[[2]]@Data))
length(ind_hvs$`Saxifraga oppositifolia`@HVList[[3]]@Data)
any(is.na(ind_hvs$`Saxifraga oppositifolia`@HVList[[3]]@Data))
length(ind_hvs$`Saxifraga oppositifolia`@HVList[[4]]@Data)
any(is.na(ind_hvs$`Saxifraga oppositifolia`@HVList[[4]]@Data))


print(head(ind_hvs$`Saxifraga oppositifolia`@HVList[[1]]@Data))
print(head(ind_hvs$`Saxifraga oppositifolia`@HVList[[2]]@Data))
print(head(ind_hvs$`Saxifraga oppositifolia`@HVList[[3]]@Data))
print(head(ind_hvs$`Saxifraga oppositifolia`@HVList[[4]]@Data))

print(ind_hvs$`Saxifraga oppositifolia`@HVList[[4]])

hv_set_list <- list()

for (i in seq_along(sp_hvs@HVList)) {
  hv_set <- hypervolume_set(sp_hvs@HVList[[i]], cavm_hv, check.memory = F)
  hv_set_list[[i]] <- hv_set
}

hv_project_results <- list()

for(i in 1:length(ind_hvs$`Saxifraga oppositifolia`@HVList)) {
  # Calculate the intersection between the species hypervolume and the region hypervolume
  hv_project <- hypervolume_project(ind_hvs$`Saxifraga oppositifolia`@HVList[[i]], biovars_cavm[[i]], type = "probability", verbose = T)
  
  # Store the result in the list
  hv_project_results[[i]] <- hv_project
}
plot_projects(hv_project_results)

plot(hv_project_results[[1]])
plot(hv_project_results[[2]])
plot(hv_project_results[[3]])
plot(hv_project_results[[4]])

length(hv_project_results[[1]]$bio_18)
values(hv_project_results[[1]])

sp_hvs@HVList[[1]]@Data[, "bio_18"]

View(sp_hvs@HVList[[1]])
View(cavm_hv)

dimnames(cavm_hv@Data)
dimnames(sp_hvs@HVList[[1]]@Data)

identical(colnames(saxi_hv@Data), names(biovars_cavm))

hv_project_results <- list()

for(i in 1:length(ind_hvs$`Saxifraga oppositifolia`@HVList)) {
  # Calculate the intersection between the species hypervolume and the region hypervolume
  hv_project <- hypervolume_project(ind_hvs$`Saxifraga oppositifolia`@HVList[[i]], biovars_cavm[[i]], type = "probability", verbose = T)
  
  # Store the result in the list
  hv_project_results[[i]] <- hv_project
}

hv_set <- hypervolume_set(ind_hvs$`Saxifraga oppositifolia`@HVList[[1]], cavm_hv@Data, check.memory=FALSE)

values(hv_set@HVList$HV1)

any(is.na(biovars_cavm[[1]]))

hv_project <- hypervolume_project(sp_hvs@HVList[[1]], biovars_cavm, type = "probability", verbose = T)

names(hv_project) <- "suitability_score"

sc <- values(hv_project[!is.na(hv_project)])

png(paste0("./outputs/data_analysis/hypervolume/species_projects/plots/hv_bio10_bio3_", names(hv_project), ".png"), width = 1920, height = 1080, pointsize = 20)

plot(hv_project, main = paste("Hypervolume_project for Saxifraga oppositifolia bio3 and 10"), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)

dev.off()
plot(hv_project)

plot(sp_hvs@HVList[[1]])
plot(cavm_hv)

hv_set <- hypervolume_set(sp_hvs@HVList[[1]], cavm_hv, check.memory = F)

hypervolume_overlap_statistics(hv_set)

# draw the Venn diagram
venn <- draw.pairwise.venn(
  area1 = length(sp_hvs@HVList[[1]]@Data),
  area2 = length(cavm_hv@Data),
  cross.area = length(intersect(sp_hvs@HVList[[1]]@Data, cavm_hv@Data)),
  category = c(sp_hvs@HVList[[1]]@Name, cavm_hv@Name),
  fill = c("blue", "red"),
  lty = "blank"
)

grid.draw(venn)

dimnames(cavm_hv@Data)
saxi <- env_data$`Saxifraga oppositifolia`
saxi_hv <- hypervolume_box(saxi,name="saxifraga oppositifolia")
hv_set <- hypervolume_set(saxi_hv, cavm_hv, check.memory=FALSE)

identical(dimnames(saxi_hv@Data), dimnames(cavm_hv@Data))




  data(penguins,package='palmerpenguins')
  penguins_no_na = as.data.frame(na.omit(penguins))
  penguins_adelie = penguins_no_na[penguins_no_na$species=="Adelie",
                                   c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
  penguins_chinstrap = penguins_no_na[penguins_no_na$species=="Chinstrap",
                                      c("bill_length_mm","bill_depth_mm","flipper_length_mm")]
  
  hv1 = hypervolume_box(penguins_adelie,name='Adelie')
  hv2 = hypervolume_box(penguins_chinstrap,name='Chinstrap')
  
  hv_set <- hypervolume_set(hv1, hv2, check.memory=FALSE)
  
  hypervolume_overlap_statistics(hv_set)

# length(unique(regions$cavm$FLOREG))
# 
# plot(regions$cavm)
# 
# regions_list <- lapply(unique(regions$cavm$BCZONE), function(zone) {
#   subset(regions$cavm, regions$cavm$BCZONE == zone)
# })
# 
# regions_list
# plot(regions_list[[3]])
# 
# test_list <- lapply(unique(regions$cavm$BCZONE), function(flora) {
#   subset(regions$cavm, regions$cavm$FLOREG == flora)
# })
# test_list
# plot(test_list[[5]])


# Get anticlockwise wkt (GBIF friendly)
# source("./hypervolume/atoms/combine_wkt_anticlockwise.R")
# anticlockwise_wkt = combine_wkt_anticlockwise(regions, max_x = T, min_x = T)
# 
# # Ask for input
# cat("Enter the bioclimatic variables you want (e.g., 'BIO4,BIO2'): \n")
# biovars_importance_sel <- biovars_importance[c(1, 3, 4, 5),]
