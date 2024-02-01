make_histogram <- function(rast, region, region.sub, region.name, plot.show = TRUE, verbose = T) {
  create_dir_if("./outputs/visualize/plots")

  
  # Figure 1A Whole CAVM
  cat(blue("Creating histogram for the entire Region\n"))
  
  if (verbose) cat("Using app() to calculate included and excluded species in each cell. \n")
  
  sp_count <- terra::app(rast, fun=function(x) c(included = sum(!is.na(x) & x == 1), excluded = sum(!is.na(x) & x == 0)) )
  
  sp_count_dt <- data.table(terra::values(sp_count, na.rm = TRUE))
  
  # Remove rows where both included and excluded are 0
  if (verbose) cat("Removing rows where both are 0. \n")
  sp_count_dt <- sp_count_dt[!(included == 0 & excluded == 0), ]
  
  if (verbose) cat("Calculating proportion of included species. \n")
  sp_count_dt$incProp <- sp_count_dt$included / (sp_count_dt$included + sp_count_dt$excluded)
  
  if (verbose) cat("Creating plot 1A. \n")
  plot1A <- ggplot(sp_count_dt, aes(x = sp_count_dt$incProp)) +
    geom_freqpoly() +
    labs(x = "Proportion of included species", y = "Number of cells", title = paste0("Potential species distribution across cells in the ", region.name)) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(plot1A)
  
  ggsave("./outputs/visualize/plots/figure-1A.png", plot1A, width = 10, height = 6)
  
  # Figure 1B different regions
  cat(blue("Creating histogram for each region in the", region.name, "\n"))
  
  if (verbose) cat("Extracting cells for the different regions. \n")
  region_cell_vals <- terra::extract(sp_count, region)
  
  region_cell_vals <- as.data.table(region_cell_vals, na.rm = TRUE)
  colnames(region_cell_vals) <- c("regionId", "included", "excluded") 
  
  region_dt <- as.data.table(region)
  
  # Add a region id to the cavm
  if (verbose) cat("Adding regionId to the", region.name, "vector. \n")
  region_dt[, regionId := .I]
  
  if (verbose) cat("Merging data tables by regionId. \n")
  sub_regions_dt <- merge(region_cell_vals, region_dt, by = "regionId")
  
  # Remove rows where both included and excluded are 0
  if (verbose) cat("Removing rows where both are 0. \n")
  sub_regions_dt <- sub_regions_dt[!(included == 0 & excluded == 0), ]
  
  if (verbose) cat("Calculating proportion of included species. \n")
  sub_regions_dt$incProp <- sub_regions_dt$included / (sub_regions_dt$included + sub_regions_dt$excluded)
  
  plot1B <- ggplot(sub_regions_dt, aes(x = sub_regions_dt$incProp, color = as.factor(sub_regions_dt[[region.sub]]))) +
    geom_freqpoly(binwidth = 0.1) +
    labs(x = "Proportion of included species", y = "Number of cells", title = paste0("Potential species distribution across cells in the floristic regions of ", region.name), color = "Floristic region") +
    scale_x_continuous(limits = c(0, NA)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(plot1B)
  
  ggsave("./outputs/visualize/plots/figure-1B.png", plot1B, width = 10, height = 6)
  
  # Plot 1C - how many species [y] have a certain proportion of cells [x]
}