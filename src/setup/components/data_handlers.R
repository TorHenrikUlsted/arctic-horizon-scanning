wrangle_dfs <- function(src, column, dynamic.name.source) {
  files <- list.files(path = src, pattern = "\\.R$")
  
  # Source all files
  lapply(paste0(src, files), source)
  
  results <- list()
  
  if (dynamic.name.source == "file") {
    for (file in files) {
      name <- sub("\\.R$", "", file)
      
      # create a name for each
      assign(name, NULL)
      
      func_name <- paste0("wrangle_", name)
      
      if (exists(func_name)) {
        assign(name, do.call(func_name, list(column = column)))
        
        # Check if the result is a data.table or data.frame
        if (!"data.table" %in% class(get(name)) && !"data.frame" %in% class(get(name))) {
          stop(paste("The result of", func_name, "is not a 'data.table' or 'data.frame'."))
        } else {
          catn(toString(class(get(name))))
        }
      }
      
      results[[name]] <- get(name)
    }
    
  } else if (dynamic.name.source == "object") {
    for (name in names(results)) {
      # Check if the object exists
      if (exists(name)) {
        # Get the object
        obj <- get(name)
        
        # Check if the object is a data.table or data.frame
        if (!"data.table" %in% class(obj) && !"data.frame" %in% class(obj)) {
          stop(paste("The object", name, "is not a 'data.table' or 'data.frame'."))
        } else {
          catn(paste("The class of the object", name, "is", class(obj)))
        }
        
        # Assign the object to the results list
        results[[name]] <- obj
      }
    }
  }
  
  return(results)
}

syncheck_dfs <- function(wrangled_dfs, column, out.dir, max.cores, verbose, counter) {
  
  synonym_lists <-  lapply(names(wrangled_dfs), function(name) {
    vebcat("Running synonym Check on", highcat(name), veb = verbose)
    
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
        vebcat("Some scientificNames are NA.", color = "nonFatalError")
        
        sp_synonyms_na <- sp_synonyms[is.na(sp_synonyms$scientificName), ]
        
        fwrite(sp_synonyms_na, paste0(out_dir, "/", child_folder, file_path, "-wfo-na.csv"), bom = T)
        
        sp_synonyms <- sp_synonyms[!is.na(sp_synonyms$scientificName), ]
        
        if (any(is.na(sp_synonyms$scientificName))) {
          vebcat("Failed at removing NA scientificNames.", color = "nonFatalError")
        } else {
          vebcat("Successfully removed NA scientificNames.", color = "proSuccess")
        } 
        
      } else {
        catn("The WFO.match result is clean.")
      }
      
      # Select best match and remove duplications
      sp_checked <- check_syn_wfo_one(sp_synonyms, paste0(out_dir, file_path))
      
      vebcat("Renaming file to resources/synonym-checked folder.", veb = verbose)
      old_name <- paste0(out_dir, "/", child_folder, "/wfo-one-uniq.csv")
      new_name <- paste0(cut_dir, "/", parent_folder, "-", child_folder, "-wfo-one.csv")
      file.rename(from = old_name, to = new_name)
      
      return(setNames(list(sp_checked), name))
      
    } else {
      catn(highcat(name), "already synonym checked.")
    }
    
  })
  
  synonym_lists <- unlist(synonym_lists, recursive = FALSE)
  
  return(synonym_lists)
}

setup_raw_data <- function(column, test = NULL, max.cores, verbose, counter) {
  vebcat("Setting up raw data.", color = "funInit")
  
  files_dir <- "./resources/synonym-checked"
  file_ext <- "-wfo-one.csv"
  files <- list.files(files_dir)
  wrangled_files <- c(
    paste0("aba-present", file_ext), 
    paste0("aba-absent", file_ext),
    paste0("ambio-present", file_ext),
    paste0("ambio-absent", file_ext),
    paste0("glonaf-species", file_ext)
  )
  
  if (!is.null(test) && length(test) > 0) {
    test <- wrangle_test(test = test, column, verbose = verbose)
    
    checklist <- test  
    
  } else {
    if (all(wrangled_files %in% files)) {
      checklist <- NULL
    } else {
      aba <- wrangle_aba(column, verbose = verbose)
      ambio <- wrangle_ambio(column, verbose = verbose)
      glonaf <- wrangle_glonaf(column, verbose = verbose)
      
      checklist <- c(aba, ambio, glonaf)
    }
  }
  
  if (!is.null(checklist)) {
    vebprint(names(checklist), verbose, "dfs added to checklist:")
    
    checked_dfs <- syncheck_dfs(
      checklist, 
      column,
      out.dir = "./outputs/setup/wrangle", 
      max.cores = max.cores, 
      verbose = verbose, 
      counter = counter
    )
    
    if (all(sapply(checked_dfs, is.null))) {
      catn("All data frames already exist.")
      
    } else {
      vebcat("Combining no-matches.", veb = verbose)
      
      combined_df <- data.frame()
      
      # Loop over each list in checked_dfs
      for(i in 1:length(checked_dfs)){
        if (!is.null(checked_dfs[[i]]$wfo_one_nomatch) && nrow(checked_dfs[[i]]$wfo_one_nomatch) > 0) {
          if (names(checked_dfs)[i] != "") {
            catn("Getting nomatches for:", highcat(names(checked_dfs)[i]))
            
            checked_dfs[[i]]$wfo_one_nomatch$dfOrigin <- names(checked_dfs)[i]
            catn("Adding an origin column with:", highcat(names(checked_dfs)[i]))
          }
          
          # rbind the wfo_one_nomatch data frame from each list
          combined_df <- rbind(combined_df, checked_dfs[[i]]$wfo_one_nomatch)
        }
      }
      
      nomatch_combined <- combined_df[!duplicated(combined_df[[column]]), ]
      
      catn("Writing combined no-matches to:", colcat("./outputs/setup/wrangle", color = "output"))
      
      if (nrow(nomatch_combined) > 0) {
        catn("There were", highcat(nrow(nomatch_combined)), "species without matches. \n")
        fwrite(nomatch_combined, "./outputs/setup/wrangle/combined-wfo-nomatch.csv", bom = T)
      } else {
        catn("There were", highcat(0), "species without any matches. \n")
      }
    }
  }
  
  vebcat("raw data setup completed successfully.", color = "funSuccess")
}