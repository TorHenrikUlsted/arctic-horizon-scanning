visualize_histogram <- function(sp_cells, region, region.sub, region.sub.color, region.name, plot.show = TRUE, verbose = T) {
  create_dir_if("./outputs/visualize/plots")
  
  # Figure 1A Whole CAVM
  cat(blue("Creating histogram for the entire Region\n"))
  
  if (verbose) cat("Creating plot 1A. \n")
  fig1A <- ggplot(sp_cells$region, aes(x = sp_cells$region$prop)) +
    geom_freqpoly() +
    labs(x = "Relative Species Distribution", y = "Cell Frequency", title = paste0("Potential species distribution in ", region.name)) +
    scale_x_continuous(limits = c(0, NA)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig1A)
  
  ggsave("./outputs/visualize/plots/figure-1A.png", fig1A, width = 10, height = 6)
  
  # Figure 1B different regions
  cat(blue("Creating histogram for each region in the", region.name, "\n"))
  
  fig1B <- ggplot(sp_cells$region, aes(x = prop, color = as.factor(sp_cells$region[[region.sub.color]]), linetype = as.factor(sp_cells$region[[region.sub]])) ) +
    geom_freqpoly(binwidth = 0.1) +
    labs(x = "Relative Species Distribution", y = "Cell Frequency", title = paste0("Potential species distribution in the floristic regions of ", region.name), color = "Country", linetype = "Floristic Province") +
    scale_x_continuous(limits = c(0, NA)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig1B)
  
  ggsave("./outputs/visualize/plots/figure-1B.png", fig1B, width = 10, height = 6)
  
  # Plot 1C - how many species [y] have a certain proportion of cells [x]
}

visualize_hotspots <- function(sp_cells, prob.value, region, region.sub, region.name, plot.show = TRUE, verbose = T) {

  inc <- sp_cells[[1]]$sp_count
  prob <- sp_cells[[2]]$sp_count

  inc[inc == 0] <- NA
  min_inc <- min(terra::values(inc$prop), na.rm = TRUE)
  max_inc <- max(terra::values(inc$prop), na.rm = TRUE)
  region_ext <- terra::ext(region)
  
  fig2 <- ggplot() +
    geom_spatvector(data = cavm, aes_string(color = paste0("as.factor(", region.sub, ")")), fill = NA) +
    scale_color_grey( guide = guide_legend(reverse = TRUE)) +
    geom_spatraster(data = inc, aes(fill = prop)) +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, breaks = c(seq(min_inc, max_inc, by = 0.2)), limits = c(min_inc, max_inc), guide = guide_legend(reverse = TRUE)) +
    labs(x = "Longitude", y = "Latitude", title = paste0("Potential species distribution in the ", region.name), fill = "Proportion", color = "Country") +
    coord_sf(xlim = c(region_ext$xmin, region_ext$xmax), ylim = c(region_ext$ymin, region_ext$ymax)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig2)
  
  ggsave("./outputs/visualize/plots/figure-2.png", plot = fig2)
  
  prob[prob == 0] <- NA
  prob_min <- min(terra::values(prob$proportion), na.rm = TRUE)
  prob_max <- max(terra::values(prob$proportion), na.rm = TRUE)
  fig3 <- ggplot() +
    geom_spatvector(data = cavm, aes_string(color = paste0("as.factor(", region.sub, ")")), fill = NA) +
    scale_color_grey( guide = guide_legend(reverse = TRUE)) +
    geom_spatraster(data = prob, aes_string(fill = prob.value)) +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, breaks = c(seq(prob_min, prob_max, by = 0.2)), limits = c(prob_min, prob_max), guide = guide_legend(reverse = TRUE)) +
    labs(x = "", y = "Latitude", title = paste0("Predicted species distribution in ", region.name), fill = "Relative Suitability Score", color = "Country") +
    coord_sf(xlim = c(region_ext$xmin, region_ext$xmax), ylim = c(region_ext$ymin, region_ext$ymax)) +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic"))
  
  if (plot.show) print(fig3)
  
  ggsave("./outputs/visualize/plots/figure-3.png", plot = fig3)
}

visualize_matrix <- function(sp_cells, region.sub, plot.show = F) {
  cells <- sp_cells$region[!is.na(prop), ]
  
  min_prop <- min(cells$prop)
  max_prop <- max(cells$prop)
  
  break_scale <- seq(min_prop, 1, by = 0.10)
  cells$prop_bin <- cut(cells$prop, breaks = break_scale, labels = FALSE)
  
  fig4 <- ggplot(cells, aes(x = prop_bin, y = as.factor(cells[[region.sub]]), fill = prop)) +
    geom_tile(color = "white") +
    scale_fill_whitebox_c(palette = "muted", na.value = NA, breaks = break_scale, guide = guide_legend(reverse = TRUE)) +
    labs(title = "Species climate suitability", x = "Species proportion", y = "Country") +
    theme_minimal() + 
    theme(plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 14, face = "bold.italic")) +
    scale_x_discrete(breaks = break_scale)
    
  if (plot.show) print(fig4)

  ggsave("./outputs/visualize/plots/figure-4.png", plot = fig4)
}

visualize_sankey <- function(df.start, df.end, df.mergeCol, plot.show = F) {
  
  print(df.start)
  print(df.end)
  print(any(df.end$included > 0))
  
  cat("Merging data frames.\n")
  merged_df <- merge(df.start, df.end, by = df.mergeCol, all = TRUE, allow.cartesian = TRUE)
  merged_df <- na.omit(merged_df)
  
  print(any(merged_df$included > 0))
  print(merged_df)
  cat("Aggregating included species with country data frames.\n")
  sum_df <- aggregate(included ~ species + country.x, data = merged_df, FUN = sum)
  
  print(sum_df)
  
  stop()
  fig5 <- ggplot(data = merged_df, aes(axis1 = country.x, axis2 = country.y, y = prop)) +
    scale_x_discrete(limits = c("included_sp country", "arctic_sub_region country"), expand = c(.1, .1)) +
    xlab("Country") +
    ylab("Prop") +
    geom_alluvium(aes(fill = species)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    theme_minimal() +
    ggtitle("Sankey plot from included_sp country to arctic_sub_region country")
  
  if (plot.show) print(fig5)
  
  ggsave("./outputs/visualize/plots/figure-5.png", plot = fig5)
  
  stop()
  
  # Figure 5 -Sankey
  ggplot(merged_inc, aes(axis1 = country, axis2 = FLOREG)) +
    geom_alluvium(aes(fill = prop)) +
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_viridis_d() +
    theme_void()
  
  ggplot(df_long, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = node)) +
    geom_sankey() +
    geom_sankey_label() +
    theme_sankey(base_size = 16)
}






