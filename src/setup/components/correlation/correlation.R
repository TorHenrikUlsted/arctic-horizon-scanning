analyze_correlation <- function(raster_stack, file.out, verbose = FALSE) {
  
  vebcat("Analyzing correlation", color = "funInit")
  
  file_out <- paste0(file.out, "/correlation-matrix.csv")
  
  if (file.exists(file_out)) {
    catn("Correlation analysis already exists at:", colcat(paste0(file.out, "/correlation-matrix.csv"), color = "output"))
    return(invisible())
  }
  
  create_dir_if(file.out)
  
  stack_corr <- terra::layerCor(raster_stack, "pearson", na.rm = T)
  
  stack_corr <- stack_corr$correlation
  
  vebprint(stack_corr, text = "Stack_corr$correlation:", veb = verbose)

  # Create correlation plot
  catn("Plotting correlation.")
  
  # Open a PNG device
  png(filename = paste0(file.out, "/corrplot-circle.png"), width = 1920, height = 1080, pointsize = 20)
  
  # Create the first plot
  corrplot(stack_corr, method = "circle", title = deparse(substitute(stack)), tl.cex = 2, cl.cex = 2, number.cex = 1, mar = c(0, 0, 2, 0))
  
  # Close the PNG device
  dev.off()
  
  # Open another PNG device
  png(filename = paste0(file.out, "/corrplot-number.png"), width = 1920, height = 1080, pointsize = 20)
  
  # Create the second plot
  corrplot(stack_corr, method = "number", title = deparse(substitute(stack)), tl.cex = 2, cl.cex = 2, number.cex = 1, mar = c(0, 0, 2, 0))
  
  # Close the PNG device
  dev.off()

  catn("Writing correlation matrix to:", colcat(file_out, color = "output"))
  
  stack_corr <- as.data.table(stack_corr)
  
  fwrite(stack_corr, file_out, row.names = T, bom = T)
  
  vebcat("Correlation analyzed successfully", color = "funSuccess")
  
  return(invisible())
}
