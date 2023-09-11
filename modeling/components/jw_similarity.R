jw_similarity_check = function(filtered_species) {
  message("------ Running Jaro-Winkler similarity check ------")
  
  ## Run a similarity check to find possible missed duplications
  simCheck_filt_sp = similarity_check(filtered_species, "scientificName", "scientificName", "jw", 0.05)
  simCheck_filt_sp = data.frame(simCheck_filt_sp)

  cat("Removing 'var.', 'x', 'susp.' from the jaro-Winkler result \n")
  ## Filter out different variances
  simCheck_filt_sp = filter_rows_after_split_text(simCheck_filt_sp, "name", "similarName", "var.")
  ## Filter out hybrids
  simCheck_filt_sp = filter_rows_around_split_text(simCheck_filt_sp, "name", "similarName", "x")
  ## Filter out different subspecies
  simCheck_filt_sp = filter_rows_after_split_text(simCheck_filt_sp, "name", "similarName", "subsp.")
  
  cat(yellow("Writing csv file to: \n", "outputs/similarity_check_outputs/jw_similarityCheck.csv \n"))
  write.csv(simCheck_filt_sp, "outputs/similarity_check_outputs/jw_similarityCheck.csv", row.names = F, fileEncoding = "UTF-8")
  
  return(simCheck_filt_sp)
}
