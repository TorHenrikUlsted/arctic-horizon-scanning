get_gbif_data = function(gbif_keys, region) {
  
  cat(cc$aquamarine("Using Keys to collect names. this may take a while \n"))
  # Set up a progress handler
  handlers("txtprogressbar")
  
  # Use with_progress instead of directly calling mclapply
  with_progress({
    # Set up a progress bar
    p = progressor(steps = length(gbif_keys))
    
    gbif_sp_data = lapply(gbif_keys, FUN = function(x) {
      # Update progress bar
      p()
      
      tryCatch({
        occ_data(
          speciesKey = x,
          hasCoordinate = T,
          geometry = region
        )$data
      }, error = function(e) {
        cat("Error with speciesKey:", x, "\n")
        cat("Error message:", e$message, "\n")
        return(NULL)
      })
    })
  })
  
  cat("Removing duplicates and merging dfs \n")
  gbif_sp_data <- lapply(gbif_sp_data, function(x) x %>% distinct(speciesKey, .keep_all = TRUE))
  gbif_sp_data <- bind_rows(gbif_sp_data)
  ## Make the species names into strings
  gbif_sp_data = data.frame(gbif_sp_data)
  
  cat("All species unique:", if(!any(duplicated(gbif_sp_data$scientificName)) == T) green("True") else red("False"), "\n")
  
  cat(green("GBIF species collected. \n"))
  
  return(gbif_sp_data)
}
