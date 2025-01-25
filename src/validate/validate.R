source_all("./src/validate/components")

validate_sequence <- function(known.file, unknown.file, out.dir = ".", vis.save.device = "jpeg", force.seq = NULL, verbose = FALSE) {
  
  vebcat("Initializing validation sequence", color = "funInit")
  
  plot_dir <- file.path(out.dir, "plots")
  
  
  res <- analyze_hypervolume_validation(
    known.file = known.file,
    unknown.file = unknown.file,
    plot.dir = plot_dir,
    vis.save.device = vis.save.device,
    verbose = verbose
  )
  
  # Add detailed analysis of no-overlap cases
  no_overlap_analysis <- analyze_no_overlap_cases(
    known.file = known.file,
    unknown.file = unknown.file,
    plot.dir = plot_dir,
    vis.save.device = vis.save.device,
    verbose = verbose
  )
  
  pattern_analysis <- analyze_patterns(
    known.file = known.file,
    unknown.file = unknown.file,
    plot.dir = plot_dir,
    vis.save.device = vis.save.device,
    verbose = verbose
  )
  
}

