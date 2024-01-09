analyze_correlation <- function(raster_stack, file.out, threshold, plot = T) {

  stack_corr <- terra::layerCor(raster_stack, "pearson", na.rm = T)
  
  stack_corr <- stack_corr$pearson
  
  create_dir_if(file.out)

  # Create correlation plot
  if (plot == T) {
    cat("Plotting correlation. \n")
    
    # Open a PNG device
    png(filename = paste0(file.out, "/corrplot_circle.png"), width = 1920, height = 1080, pointsize = 20)
    
    # Create the first plot
    corrplot(stack_corr, method = "circle", title = deparse(substitute(stack)), tl.cex = 2, cl.cex = 2, number.cex = 1, mar = c(0, 0, 2, 0))
    
    # Close the PNG device
    dev.off()
    
    # Open another PNG device
    png(filename = paste0(file.out, "/corrplot_number.png"), width = 1920, height = 1080, pointsize = 20)
    
    # Create the second plot
    corrplot(stack_corr, method = "number", title = deparse(substitute(stack)), tl.cex = 2, cl.cex = 2, number.cex = 1, mar = c(0, 0, 2, 0))
    
    # Close the PNG device
    dev.off()
    
  } else {
    cat("Correlation plotting skipped. \n")
  }
  
  fwrite(stack_corr, paste0(file.out, "/correlation_matrix.csv"), row.names = T, bom = T)
  
  return(stack_corr)
}