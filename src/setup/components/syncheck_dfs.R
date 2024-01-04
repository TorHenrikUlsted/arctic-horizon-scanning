syncheck_dfs <- function(wrangled_dfs, column, out.dir, max.cores, verbose, counter) {

  synonym_lists <-  lapply(names(wrangled_dfs), function(name) {
    if (verbose) cat("Running synonym Check on", cc$lightSteelBlue(name), "\n")
    
    split_name <- strsplit(name, "_")[[1]]
    parent_folder <- split_name[1]
    child_folder <- split_name[2]
    
    out_dir <- paste0(out.dir, "/", parent_folder)
    create_dir_if(out_dir)
    
    cut_dir <- paste0("./resources/synonym-checked")
    create_dir_if(cut_dir)
    
    file_path <- paste0("/", child_folder)
    
    if(!file.exists(paste0(cut_dir, "/", parent_folder, "-", child_folder, "-wfo-one.csv"))) {
      
      # Run synonym check on the species
      sp_synonyms <- check_syn_wfo(
        checklist = wrangled_dfs[[name]],
        column = column,
        folder = paste0(out_dir, "/", child_folder),
        max.cores = max.cores,
        verbose = verbose,
        counter = counter
        )
      
      
      if (any(is.na(sp_synonyms$scientificName))) {
        cat(cc$lightCoral("Some scientificNames are NA. \n"))
        
        sp_synonyms_na <- sp_synonyms[is.na(sp_synonyms$scientificName), ]
        
        fwrite(sp_synonyms_na, paste0(out_dir, "/", child_folder, file_path, "-wfo-na.csv"), bom = T)
        
        sp_synonyms <- sp_synonyms[!is.na(sp_synonyms$scientificName), ]
        
        if (any(is.na(sp_synonyms$scientificName))) cat(cc$lightCoral("Failed at removing NA scientificNames. \n")) else cat(cc$lightGreen("Successfully removed NA scientificNames. \n"))
      } else {
        cat("The WFO.match result is clean. \n")
      }
      
      # Select best match and remove duplications
      sp_checked <- check_syn_wfo_one(sp_synonyms, paste0(out_dir, file_path))
      
      if (verbose) cat("Renaming file to resources/synonym-checked folder. \n")
      old_name <- paste0(out_dir, "/", child_folder, "-wfo-one.csv")
      new_name <- paste0(cut_dir, "/", parent_folder, "-", child_folder, "-wfo-one.csv")
      file.rename(from = old_name, to = new_name)
      
      return(setNames(list(sp_checked), name))
      
    } else {
      cat(cc$lightSteelBlue(name), "already synonym checked. \n")
    }

  })
  
  synonym_lists <- unlist(synonym_lists, recursive = FALSE)
  
  return(synonym_lists)
}
