plot_hypervolumes <- function(hv_list, out.dir) {
  vebcat("Plotting hypervolumes.", color = "funInit")
  
  for (i in seq_along(hv_list)) {
    hv <- hv_list[[i]]
    
    create_dir_if(out.dir)
    
    catn("Plotting", highcat(names(hv_list[[i]])))
    
    png(paste0(out.dir, "/hv_", names(hv_list)[i], ".png"), width = 800, height = 800, pointsize = 20)
    
    plot(hv, main = paste("Hypervolume for species", names(hv_list)[i]), cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
    
    dev.off()
  }
  
  vebcat("Plotting hypervolumes completed successfully", color = "funSuccess")
}