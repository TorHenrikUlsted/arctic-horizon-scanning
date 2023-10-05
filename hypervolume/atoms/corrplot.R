plot_correlation <- function(list) {
  # Check dimensions of each raster object
  dims <- lapply(list, function(x) dim(x))

  # Check if all dimensions are the same
  cat("All matrix dimensions are equal: ", if (length(unique(dims)) == 1) green("True") else red("false"), "\n")
  cat("Matrix sample: ", as.character(dim(list[[1]])), "\n")

  # Load list into dataframe
  biovars_df <- data.frame(lapply(list, function(x) as.vector(values(x))))
  names(biovars_df) <- paste0("Var", seq_along(list))
  # Remove NA values
  biovars_df <- na.omit(biovars_df)
  cat("All NAs removed: ", if (any(is.na(biovars_df))) red("False") else green("True"), "\n")

  # Calculate the correlation matrix
  cor_mat <- cor(biovars_df, use = "pairwise.complete.obs")

  # Create correlation plot
  corrplot(cor_mat, method = "circle")
}