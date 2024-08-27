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
      if (exists(name)) {
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

  if (grepl("test", fun.name)) {
    name <- sub("^wrangle_test_", "", fun.name)
    dir <- paste0("./outputs/setup/wrangle/test/", name)
  } else {
    dir <- paste0("./outputs/setup/wrangle/", name)
    present_out <- paste0(dir, "/", name, "-present.csv")
    absent_out <- paste0(dir, "/", name, "-absent.csv")
  }
  
  vebcat("List name:", highcat(name), veb = verbose)
  vebcat("Directory:", highcat(dir), veb = verbose)
  
  create_dir_if(dir)
  present_out <- paste0(dir, "/", name, "-present.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
  print(present_out)
  print(absent_out)
  
  if (!file.exists(present_out)) present_out <- NULL
  if (!file.exists(absent_out)) absent_out <- NULL
  
  if ((!is.null(present_out) || !is.null(absent_out))) {
    catn(highcat(name), "already wrangled, loading files..")
    
    print("")

    res <- list()

    if (!is.null(present_out)) res$present <- fread(present_out, sep = "\t")
    if (!is.null(absent_out)) res$absent <- fread(absent_out, sep = "\t")

    catn("files loaded.")
  } else {
    vebcat("Initiating", name, "wrangling protocol.", color = "funInit")

    res <- fun(name = name, column = column, verbose = verbose)
    
    vebcat(name, "wrangling protocol successfully completed.", color = "funSuccess")
  }

  name <- sub("^wrangle_", "", fun.name)
  
  result <- list()
  if (!is.null(res$present)) result[[paste0(name, "_present")]] <- res$present
  if (!is.null(res$absent)) result[[paste0(name, "_absent")]] <- res$absent
  
  return(result)
}

wrangle_all <- function(column, verbose = FALSE) {
  tryCatch(
    {
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
    },
    error = function(e) {
      vebcat("Error occurred when trying to wrangle all file.", color = "fatalError")
      stop(e$message)
    }
  )

  return(results)
}

syncheck_dfs <- function(wrangled_dfs, column, out.dir, cores.max, verbose, counter) {
  synonym_lists <- lapply(names(wrangled_dfs), function(name) {
    vebcat("Running synonym Check on", highcat(name), veb = verbose)

    split_name <- strsplit(name, "_")[[1]]
    parent_folder <- split_name[1]
    child_folder <- split_name[2] # present / absent
    
    out_dir <- paste0(out.dir, "/", parent_folder)
    
    out_dir_child <- paste0(out_dir, "/", child_folder)
    
    create_dir_if(out_dir_child)

    check_file <- paste0(out_dir_child, "/wfo-one-clean.csv")

    file_name <- paste0(parent_folder, "-", child_folder, "-wfo")
    
    if (!file.exists(check_file)) {
      # Run synonym check on the species
      sp_synonyms <- check_syn_wfo(
        checklist = wrangled_dfs[[name]],
        column = column,
        out.dir = out_dir_child,
        cores.max = cores.max,
        verbose = verbose,
        counter = counter
      ) # returns list of clean and mismatched info
      
      # Select best match
      sp_checked <- check_syn_wfo_one(
        wfo.match.dt = sp_synonyms$clean,
        column = column,
        out.dir = out_dir_child,
        verbose = verbose
      )
      
      vebprint(sp_checked, text = "WFO.one checked list:")
      
      manual_checks <- sum(
        nrow(sp_synonyms$mismatch), 
        nrow(sp_checked$nomatch), 
        nrow(sp_checked$na), 
        na.rm = TRUE
      )
      
      lost_diff <- (nrow(wrangled_dfs[[name]]) - nrow(sp_checked$clean))
      
      md_dt <- data.table(
        wrangle = nrow(wrangled_dfs[[name]]), 
        match = nrow(wfo_match_dt),
        one = nrow(sp_checked$raw),
        result = nrow(sp_checked$clean),
        lost = lost_diff,
        manual = manual_checks,
        mismatch = nrow(sp_synonyms$mismatch),
        nomatch = nrow(sp_checked$nomatch),
        na = nrow(sp_checked$na),
        duplicate = nrow(sp_checked$duplicate)
      )
      
      mdwrite(
        config$files$post_seq_md,
        text = paste0("3;Standardization results ", name, " :"),
        data = md_dt
      )
      
      sp_checked$raw = NULL
      
      # Choose approach
      if (config$simulation$approach == "precautionary") {
        catn("Using precautionary method to remove infraSpecificEpithets")
        orig_n <- nrow(sp_checked$clean)
        sp_checked$clean[, scientificName := NULL]
        sp_checked$clean[, scientificName := fifelse(is.na(genus) | is.na(specificEpithet), 
                                       NA_character_, 
                                       paste(trimws(genus), " ", trimws(specificEpithet)))]
        
        sp_checked$clean <- unique(sp_checked$clean, by = "scientificName")
        new_n <- nrow(sp_checked$clean)
        
        catn(highcat(orig_n - new_n), "infraspecificEpithets removed")
        
        mdwrite(
          config$files$post_seq_md,
          text = paste0(
            "3;Standardization precautionary conversion\n\n",
            "Removed ", orig_n - new_n, " infraspecificEpithets"
          )
        )
        
        if (file.exists(check_file)) file.remove(check_file)
        fwrite(sp_checked$clean, check_file, bom = TRUE)
      }
      
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

  files_dir <- "./outputs/setup/wrangle"
  
  create_dir_if(files_dir)
  
  checklist <- wrangle_all(column = column, verbose = verbose)
  
  if (!is.null(checklist)) {
    vebprint(names(checklist), verbose, "dfs added to checklist:")

    checked_dts <- syncheck_dfs(
      checklist,
      column,
      out.dir = files_dir,
      cores.max = cores.max,
      verbose = verbose,
      counter = counter
    )
    
    if (all(sapply(checked_dts, is.null))) {
      catn("All data frames already exist.")
    } else {
      
      vebcat("Combining no-matches.", veb = verbose)

      combined_dt <- data.table()

      # Loop over each list in checked_dfs
      for (dt_name in names(checked_dts)) {
        wfo_one <- checked_dts[[dt_name]]
        if (!is.null(wfo_one$clean)) wfo_one$clean <- NULL
        if (!is.null(wfo_one$clean)) wfo_one$duplicate <- NULL
        
        vebprint(wfo_one, verbose, paste0("Synonym checked data table ", dt_name, " :"))
        
        dts_to_combine <- list()
        
        for (component in c("mismatch", "nomatch", "na")) {
          if (!is.null(wfo_one[[component]]) && nrow(wfo_one[[component]]) > 0) {
            catn("Appending", component, "for:", highcat(dt_name))
            wfo_one[[component]][, listOrigin := dt_name]
            wfo_one[[component]][, resultType := component]
            dts_to_combine <- c(dts_to_combine, list(wfo_one[[component]]))
          }
        }
        
        if (length(dts_to_combine) > 0) {
          combined_dt <- rbindlist(list(combined_dt, dts_to_combine), fill = TRUE)
        }
        
      }

      manual_combined <- unique(combined_dt, by = column)

      catn("Writing combined no-matches to:", colcat(files_dir, color = "output"))
      
      mdwrite(
        config$files$post_seq_md,
        text = paste0("3;Combined nomatches: **", nrow(nomatch_combined), "**")
      )

      if (nrow(nomatch_combined) > 0) {
        catn("There were", highcat(nrow(nomatch_combined)), "species without matches. \n")
        vebcat("These need to be checked manually.", color = "indicator")
        fwrite(nomatch_combined, paste0(files_dir, "/combined-wfo-nomatch.csv"), bom = T)
      } else {
        catn("There were", highcat(0), "species without any matches. \n")
      }
    }
  }

  vebcat("raw data setup completed successfully.", color = "funSuccess")
}
