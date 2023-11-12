get_correlation <- function(list) {
  # Check dimensions of each raster object
  dims <- lapply(list, function(x) dim(x))

  # Check if all dimensions are the same
  cat("All matrix dimensions are equal: ", if (length(unique(dims)) == 1) green("True") else red("false"), "\n")
  cat("Matrix sample: ", as.character(dim(list[[1]])), "\n")
  
  # cat("Adding biovars descriptive names \n")
  # desc_names <- c("AnnualMeanTemp", "MeanDiurnalRange", "Isothermality", "TempSeasonality", "MaxTempWarmestMonth", "MaxTempColdestMonth", "TempAnnualRange", "MeanTempWettestQuarter", "MeanTempDriestQuarter", "MeanTempWarmestQuarter", "MeanTempColdestQuarter", "AnnualPrecipitation", "PrecipitationWettestMonth", "PrecipitationDriestMonth", "PrecipitationSeasonality", "PrecipitationWettestQuarter", "PrecipitationDriestQuarter", "PrecipitationWarmestQuarter", "PrecipitationColdestQuarter")

  # Load list into dataframe
  biovars_df <- data.frame(lapply(list, function(x) as.vector(values(x))))
  
  names(biovars_df) <- paste0("BIO",seq_along(list))
  #names(biovars_df) <- desc_names
  
  # Remove NA values
  biovars_df <- na.omit(biovars_df)
  
  cat("All NAs removed: ", if (any(is.na(biovars_df))) red("False") else green("True"), "\n")

  # Calculate the correlation matrix
  cor_mat <- cor(biovars_df, use = "pairwise.complete.obs")

  # Create correlation plot
  corrplot(cor_mat, method = "circle")
  
  # Calculate the absolute correlations with all other variables
  abs_correlations <- abs(cor_mat)
  
  # Get the mean of absolute correlations for each variable
  mean_abs_correlations <- rowMeans(abs_correlations)
  
  # Order the variables by mean of absolute correlations in descending order
  important_biovars <- order(mean_abs_correlations, decreasing = TRUE)
  
  # Make data frame
  important_biovars_df <- data.frame(Variable = names(mean_abs_correlations)[important_biovars], Importance = mean_abs_correlations[important_biovars])
  
  # Print the names of the variables in order of importance
  cat("Most important BioVariables in order:", names(mean_abs_correlations)[important_biovars], "\n")
  
  fdir <- "./outputs/hypervolume/correlation/"
  
  if (!dir.exists(fdir)) dir.create(fdir, recursive = T)
  
  fwrite(important_biovars_df, paste0(fdir, "most_important_vars.txt"), sep="\t", row.names = F, bom = T)
  
  return(important_biovars_df)
}