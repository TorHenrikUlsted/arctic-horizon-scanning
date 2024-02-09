visualize_histogram <- function(rast, region, region.sub, region.name, plot.show = TRUE, verbose = T) {
  create_dir_if("./outputs/visualize/plots")

  
  # Figure 1A Whole CAVM
  cat(blue("Creating histogram for the entire Region\n"))
  
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