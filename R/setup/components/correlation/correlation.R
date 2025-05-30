analyze_correlation <- function(raster_stack, dir.out, verbose = FALSE) {
  cli::cli_h1("Correlation Analysis")

  corr_file <- file.path(dir.out, "correlation-matrix.csv")
  corr_p_file <- file.path(dir.out, "correlation-pvalues.csv")
  pca_file <- file.path(dir.out, "correlation-variance-explained.csv")

  if (file.exists(corr_file) && file.exists(corr_p_file) && file.exists(pca_file)) {
    cli::cli_alert_warning("Correlation analysis already exists at: {.file {dir.out}}")
    return(invisible())
  }

  cli::cli_progress_step("Creating output directory")
  create_dir_if(dir.out)

  cli::cli_progress_step("Computing correlations")
  stack_corr <- terra::layerCor(raster_stack, "pearson", na.rm = T)

  corr_matrix <- stack_corr$correlation

  # Add significance testing
  cli::cli_progress_step("Calculating statistical significance")
  nobs <- stack_corr$n[2, 1] # Take the second number because the first is NAN all numbers are the same
  t_stat <- corr_matrix * sqrt((nobs - 2) / (1 - corr_matrix^2))
  p_values <- 2 * pt(-abs(t_stat), df = nobs - 2)

  # Add PCA
  cli::cli_progress_step("Performing PCA analysis")
  # Extract values and remove NAs before PCA
  vals <- values(raster_stack)
  vals <- na.omit(vals)
  pca_result <- prcomp(vals, scale. = TRUE)
  var_explained <- summary(pca_result)$importance[2, ]
  cumulative_var <- summary(pca_result)$importance[3, ]

  # Plot creation
  cli::cli_progress_step("Creating correlation plots")

  cli::cli_alert("Saving circle plot")
  png(
    filename = file.path(dir.out, "corrplot-circle.png"),
    width = 1920, height = 1080, pointsize = 20
  )
  corrplot(corr_matrix,
    method = "circle",
    title = deparse(substitute(stack)),
    tl.cex = 2, cl.cex = 2, number.cex = 1,
    mar = c(0, 0, 2, 0)
  )
  dev.off()

  cli::cli_alert("Saving number plot")
  png(
    filename = file.path(dir.out, "corrplot-number.png"),
    width = 1920, height = 1080, pointsize = 20
  )
  corrplot(corr_matrix,
    method = "number",
    title = deparse(substitute(stack)),
    tl.cex = 2, cl.cex = 2, number.cex = 1,
    mar = c(0, 0, 2, 0)
  )
  dev.off()

  # Plot p-values
  cli::cli_alert("Saving p-values plot")
  png(
    filename = file.path(dir.out, "p-values-plot.png"),
    width = 1920, height = 1080, pointsize = 20
  )
  corrplot(p_values,
    method = "color",
    is.corr = FALSE,
    title = "P-Values",
    tl.cex = 2, cl.cex = 2, number.cex = 1,
    mar = c(0, 0, 2, 0),
    insig = "blank"
  )
  dev.off()

  # Plot PCA variance explained
  cli::cli_alert("Saving PCA variance explained plot")
  png(
    filename = file.path(dir.out, "pca-variance-explained-plot.png"),
    width = 1920, height = 1080, pointsize = 20
  )
  barplot(var_explained,
    main = "PCA Variance Explained",
    xlab = "Principal Components", ylab = "Variance Explained",
    col = "blue"
  )
  dev.off()

  # Scree plot
  png(filename = file.path(dir.out, "scree_plot.png"), width = 1920, height = 1080, pointsize = 20)
  plot(pca_result, type = "lines", main = "Scree Plot")
  dev.off()

  # Cumulative variance plot
  png(filename = file.path(dir.out, "cumulative_variance_plot.png"), width = 1920, height = 1080, pointsize = 20)
  plot(cumulative_var, type = "b", main = "Cumulative Variance Explained", xlab = "Principal Components", ylab = "Cumulative Variance Explained")
  dev.off()

  # Biplot
  png(filename = file.path(dir.out, "biplot.png"), width = 1920, height = 1080, pointsize = 20)
  biplot(pca_result, main = "PCA Biplot")
  dev.off()

  # Loadings plot for the first two principal components
  loadings <- pca_result$rotation[, 1:4]
  png(filename = file.path(dir.out, "loadings_plot.png"), width = 1920, height = 1080, pointsize = 20)
  barplot(t(loadings), beside = TRUE, legend = colnames(loadings), main = "Loadings Plot", xlab = "Variables", ylab = "Loadings")
  dev.off()

  # Save results
  cli::cli_progress_step("Writing results to files")

  corr_matrix_dt <- as.data.table(corr_matrix)
  fwrite(corr_matrix_dt, corr_file, row.names = T, bom = T)

  p_values_dt <- as.data.table(p_values)
  fwrite(p_values_dt, corr_p_file,
    row.names = T, bom = T
  )

  var_explained_dt <- data.table(
    PC = seq_along(var_explained),
    Variance_Explained = var_explained,
    Cumulative_Variance = cumulative_var
  )

  # Sort by Variance_Explained in descending order
  sorted_var_explained_dt <- var_explained_dt[order(-Variance_Explained)]
  fwrite(sorted_var_explained_dt, file.path(dir.out, "sorted_variance_explained.csv"))

  # Print summary statistics
  # Print summary statistics
  cli::cli_h2("Summary Statistics")
  cli::cli_alert_success("Number of observations: {.val {nobs}}")
  cli::cli_alert_success("Total variance explained by PCs: {.val {round(sum(var_explained)*100, 2)}}%")
  cli::cli_alert_success("First 4 PCs explain: {.val {round(sum(var_explained[1:4])*100, 2)}}%")

  # Detailed variance explained by each of the top 4 PCs
  for (i in 1:4) {
    cli::cli_alert_info("PC{i}: {.val {round(var_explained[i]*100, 2)}}% variance explained")
  }

  # Loadings of the top 4 principal components
  cli::cli_h3("Top 4 Principal Components Loadings")
  top_pcs <- 4
  loadings <- pca_result$rotation[, 1:top_pcs]
  bio_names <- names(raster_stack)
  summary_loadings <- apply(loadings, 1, function(x) sum(abs(x)))

  # Combine biovariable names with their summarized loadings
  biovariable_summary <- data.table(
    Biovariable = bio_names,
    Summarized_Loadings = summary_loadings
  )

  # Sort biovariables by summarized loadings
  sorted_biovariable_summary <- biovariable_summary[order(-Summarized_Loadings)]
  fwrite(sorted_biovariable_summary, file.path(dir.out, "sorted_biovariables_loadings.csv"), bom = TRUE)

  cli::cli_h3("Summarized Loadings for Top 4 Principal Components")
  print(as.data.table(sorted_biovariable_summary))

  cli::cli_alert_success("Analysis completed successfully!")

  return(invisible(list(
    correlation = corr_matrix,
    p_values = p_values,
    variance_explained = var_explained
  )))
}
