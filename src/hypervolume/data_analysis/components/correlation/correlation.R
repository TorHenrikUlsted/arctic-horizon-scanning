analyze_correlation <- function(raster_stack, file.out, threshold, plot = T, verbose = F) {
  
  vebcat("Analyzing correlation", color = "funInit")
  
  if (file.exists(paste0(file.out, "/correlation_matrix.csv"))) {
    catn("Correlation analysis already exists at:",  paste0(file.out, "/correlation_matrix.csv"))
    return()
  }
  
  create_dir_if(file.out)
  
  stack_corr <- terra::layerCor(raster_stack, "pearson", na.rm = T)
  
  stack_corr <- stack_corr$correlation
  
  vebprint(stack_corr, text = "Stack_corr$correlation:")

  # Create correlation plot
  if (plot == T) {
    vebcat("Plotting correlation.", veb = verbose)
    
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
    vebcat("Correlation plotting skipped.", veb = verbose)
  }
  
  file_out <- paste0(file.out, "/correlation_matrix.csv")
  
  catn("Writing correlation matrix to:", colcat(file_out, color = "output"))
  
  fwrite(stack_corr, file_out, row.names = T, bom = T)
  
  vebcat("Correlation analyzed successfully", color = "funSuccess")
  
  return(stack_corr)
}