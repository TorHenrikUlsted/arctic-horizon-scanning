source_all("./R/validate/components")

validate_sequence <- function(known.dir, unknown.dir, out.dir = ".", vis.save.device = "jpeg", force.seq = NULL, verbose = FALSE) {
  vebcat("Initializing validation sequence", color = "funInit")

  plot_dir <- file.path(out.dir, "plots")
  known_file <- file.path(known.dir, "stats.csv")
  unknown_file <- file.path(unknown.dir, "stats.csv")

  res <- analyze_hypervolume_validation(
    known.dir = known.dir,
    unknown.dir = unknown.dir,
    example = config$simulation$example,
    plot.dir = plot_dir,
    vis.save.device = vis.save.device,
    verbose = verbose
  )

  # vebprint(res, text = "Res:")

  # Add detailed analysis of no-overlap cases
  no_overlap_analysis <- analyze_no_overlap_cases(
    known.file = known_file,
    unknown.file = unknown_file,
    plot.dir = plot_dir,
    vis.save.device = vis.save.device,
    verbose = verbose
  )

  pattern_analysis <- analyze_patterns(
    known.file = known_file,
    unknown.file = unknown_file,
    permutations = 999,
    plot.dir = plot_dir,
    vis.save.device = vis.save.device,
    verbose = verbose
  )


  # vebprint(no_overlap_analysis, text = "No overlap Analysis:")
  # vebprint(pattern_analysis, text = "Pattern Analysis:")
}
