check_wrangled_dfs <- function(src, column) {
  # Wrangle all lists
  wrangled_lists <- wrangle_dfs(src, column)

  # Synonym Check all the lists separately.
  synonym_lists <- lapply(names(wrangled_lists), function(name) {
    split_name <- strsplit(name, "_")[[1]]
    directory <- "setup"
    parent_folder <- split_name[1]
    child_folder <- split_name[2]

    # Run synonym check on the species
    sp_synonyms <- check_syn_wfo(wrangled_lists[[name]], column, paste(directory, parent_folder, child_folder, sep = "/"))

    # Select best match and remove duplications
    sp_checked <- check_syn_wfo_one(sp_synonyms, column, paste(directory, parent_folder, child_folder, sep = "/"))

    return(sp_checked)
  })

  return(synonym_lists)
}
