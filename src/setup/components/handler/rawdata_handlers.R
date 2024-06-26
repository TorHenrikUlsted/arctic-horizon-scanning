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

wrangle_if <- function(fun.name, column, verbose = FALSE) {
  fun <- get(fun.name)
  
  name <- sub("^wrangle_", "", fun.name)
  
  vebcat("List name:", highcat(name), veb = verbose)
  
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)
  
  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  
  present_out <- paste0(dir, "/", name, "-present.csv")
  if (!file.exists(present_out)) present_out = NULL
  
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  if (!file.exists(absent_out)) absent_out = NULL
  
  if ((!is.null(present_out) || !is.null(absent_out))) {
    catn(highcat(name), "already wrangled, loading files..")
    
    res <- list()
    
    if (!is.null(present_out)) res$present <- fread(present_out, sep = "\t")
    if (!is.null(absent_out)) res$absent <- fread(absent_out, sep = "\t")
    
    catn("files loaded.")
  } else {
    vebcat("Initiating", name, "wrangling protocol.", color = "funInit")
    
    res <- fun(name = name, column = column, verbose = verbose)
    
    vebcat(name, "wrangling protocol successfully completed.", color = "funSuccess")
  }
  
  result <- list()
  if (!is.null(res$present)) result[[paste0(name, "_present")]] <- res$present
  if (!is.null(res$absent)) result[[paste0(name, "_absent")]] <- res$absent
    
  return(result)
}

wrangle_all <- function(column, verbose = FALSE) {
  tryCatch({
    all_funs <- ls(envir = .GlobalEnv)[sapply(mget(ls(envir = .GlobalEnv), envir = .GlobalEnv), is.function)]
    
    wrangle_funs <- grep("^wrangle_", all_funs, value = TRUE)
    
    wrangle_funs <- setdiff(wrangle_funs, c("wrangle_all", "wrangle_dfs", "wrangle_template", "wrangle_if"))
    
    vebprint(wrangle_funs, text = "all Wrangle functions:")
    
    results <- list()
    
    for (fun_name in wrangle_funs) {
      name <- sub("^wrangle_", "", fun_name)
      
      vebcat("Wrangling function:", fun_name, veb = verbose)
      
      res <- wrangle_if(fun.name = fun_name, column = column, verbose = verbose)
      
      for (sublist_name in names(res)) {
        results[[paste0(name, "$", sublist_name)]] <- res[[sublist_name]]
      }
      
      names(results) <- gsub(".*\\$", "", names(results))
    }
  }, error = function(e) {
    vebcat("Error occurred when trying to wrangle all file.", color = "fatalError")
    stop(e$message)
  })
  
  return(results)
}

syncheck_dfs <- function(wrangled_dfs, column, out.dir, cores.max, verbose, counter) {
  
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
        cores.max = cores.max,
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

setup_raw_data <- function(column, cores.max = 1, verbose = FALSE, counter = 1) {
  vebcat("Setting up raw data.", color = "funInit")
  
  files_dir <- "./resources/synonym-checked"
  file_ext <- "-wfo-one.csv"
  files <- list.files(files_dir)
  
  checklist <- wrangle_all(column = column, verbose = verbose)
  
  if (!is.null(checklist)) {
    vebprint(names(checklist), verbose, "dfs added to checklist:")
    
    checked_dfs <- syncheck_dfs(
      checklist, 
      column,
      out.dir = "./outputs/setup/wrangle", 
      cores.max = cores.max, 
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
        vebcat("These need to be checked manually.", color = "indicator")
        fwrite(nomatch_combined, "./outputs/setup/wrangle/combined-wfo-nomatch.csv", bom = T)
      } else {
        catn("There were", highcat(0), "species without any matches. \n")
      }
    }
  }
  
  vebcat("raw data setup completed successfully.", color = "funSuccess")
}
