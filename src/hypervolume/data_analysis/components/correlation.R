get_correlation <- function(stack, threshold, plot = T) {
  # Check dimensions of each raster object
  dims <- dim(stack)
  
  # Get the dimensions of each layer in the stack
  layer_dims <- lapply(1:nlyr(stack), function(i) dim(stack[[i]]))
  
  # Check if all layers have the same number of rows and columns
  all_dims_equal <- length(unique(layer_dims)) == 1
  
  cat("All layers have the same dimensions: ", if (all_dims_equal) green("TRUE") else red("FALSE"), "\n")
  cat("Matrix sample: ", as.character(dim(stack)), "\n")

  # Load list into dataframe
  biovars_df <- as.data.frame(stack, xy=FALSE)
  
  names(biovars_df) <- paste0("bio_", gsub(".*_", "", names(stack)))
  
  # Remove NA values
  biovars_df <- na.omit(biovars_df)
  
  cat("All NAs removed: ", if (any(is.na(biovars_df))) red("False") else green("True"), "\n")

  # Calculate the correlation matrix
  cor_mat <- cor(biovars_df, use = "pairwise.complete.obs")
  
  print(cor_mat)

  # Create correlation plot
  if (plot == T) {
    cat("Plotting correlation. \n")
    
    # Open a PNG device
    png(filename = "./outputs/data_analysis/correlation/corrplot_circle.png", width = 1920, height = 1080, pointsize = 20)
    
    # Create the first plot
    corrplot(cor_mat, method = "circle", title = deparse(substitute(stack)), tl.cex = 2, cl.cex = 2, number.cex = 1, mar = c(0, 0, 2, 0))
    
    # Close the PNG device
    dev.off()
    
    # Open another PNG device
    png(filename = "./outputs/data_analysis/correlation/corrplot_number.png", width = 1920, height = 1080, pointsize = 20)
    
    # Create the second plot
    corrplot(cor_mat, method = "number", title = deparse(substitute(stack)), tl.cex = 2, cl.cex = 2, number.cex = 1, mar = c(0, 0, 2, 0))
    
    # Close the PNG device
    dev.off()
    
  } else {
    cat("Correlation plotting skipped. \n")
  }
  
  fdir <- "./outputs/data_analysis/correlation/"
  
  if (!dir.exists(fdir)) dir.create(fdir, recursive = T)
  
  fwrite(cor_mat, paste0(fdir, "correlation_matrix.csv"), row.names = T, bom = T)
  
  return(cor_mat)
}