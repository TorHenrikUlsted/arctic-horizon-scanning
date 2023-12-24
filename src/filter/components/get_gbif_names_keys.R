get_gbif_names_keys <- function(region, counts, cores.max = 1, verbose = F) {
  cat(blue("Initiating GBIF name by key acquisition protocol. \n"))

    handlers("txtprogressbar")
  
  # Use with_progress instead of directly calling mclapply
  with_progress({
    # Set up a progress bar
    p <- progressor(steps = length(counts))
    
    gbif_names <- mclapply(counts, FUN = function(x) {
      # Update progress bar
      p()
      
      occ_count(
        speciesKey = x,
        facet = "scientificName",
        hasCoordinate = T
      )$scientificName
    }, mc.cores = cores.max)
  })
  
  ## Make the species names into strings
  gbif_names <- as.character(do.call(c, gbif_names))
  gbif_names <- data.frame(gbif_names = trimws(gbif_names))
  
  names_savefile <- paste0("./outputs/filter/gbif-acquisition/gbif-names-", region, ".csv")
  cat(yellow("Writing out gbif_names to: \n", names_savefile, "\n"))
  gbif_names <- set_df_utf8(gbif_names)
  fwrite(gbif_names, names_savefile, row.names = F, bom = T)
  
  cat(cc$lightGreen("GBIF name by key acquisition protocol completed successfully. \n"))
  return(gbif_names)
}