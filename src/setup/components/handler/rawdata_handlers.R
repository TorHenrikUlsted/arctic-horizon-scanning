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
  
  if (!file.exists(present_out)) present_out <- NULL
  
  if (!file.exists(absent_out)) absent_out <- NULL

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

    check_file <- paste0(out_dir_child, "/wfo-one-uniq.csv")

    file_name <- paste0(parent_folder, "-", child_folder, "-wfo")
    
    if (!file.exists(check_file)) {
      # Run synonym check on the species
      sp_synonyms <- check_syn_wfo(
        checklist = wrangled_dfs[[name]],
        column = column,
        folder = out_dir_child,
        cores.max = cores.max,
        verbose = verbose,
        counter = counter
      )

      if (any(is.na(sp_synonyms$scientificName))) {
        vebcat("Some scientificNames are NA.", color = "nonFatalError")

        sp_synonyms_na <- sp_synonyms[is.na(sp_synonyms$scientificName), ]

        fwrite(sp_synonyms_na, paste0(out_dir_child, file_name, "-na.csv"), bom = T)

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
      sp_checked <- check_syn_wfo_one(
        sp_synonyms, 
        folder = out_dir_child
      )
      
      origin_dt <- length(unique(wrangled_dfs[[name]], by = column)[[column]])
      wfo_match_dt <- length(unique(rbindlist(sp_synonyms), by = "scientificName")$scientificName)
      wfo_one_dt <- length(unique(rbindlist(sp_checked), by = "scientificName")$scientificName)
      
      md_dt <- data.table(
        wrangledSpecies = origin_dt, 
        wfoMatchSpecies = wfo_match_dt,
        wfoOneSpecies = wfo_one_dt
      )
      
      mdwrite(
        config$files$post_seq_md,
        text = paste0("2;WFO.one result species number for ", name, " :"),
        data = md_dt
      )
      
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

    checked_dfs <- syncheck_dfs(
      checklist,
      column,
      out.dir = files_dir,
      cores.max = cores.max,
      verbose = verbose,
      counter = counter
    )
    
    files <- get_files(files_dir, include.files = "wfo-one-uniq.csv")
    print(files)
    
    # Add number of species per df
    

    if (all(sapply(checked_dfs, is.null))) {
      catn("All data frames already exist.")
    } else {
      vebcat("Combining no-matches.", veb = verbose)

      combined_df <- data.frame()

      # Loop over each list in checked_dfs
      for (i in 1:length(checked_dfs)) {
        dt <- checked_dfs[[i]]
        dt_name <- names(checked_dfs)[i]
        vebprint(dt, verbose, paste0("Synonym checked data table ", dt_name, " :"))
        
        if (!is.null(dt$wfo_one_nomatch) && nrow(dt$wfo_one_nomatch) > 0) {
          if (dt_name != "") {
            catn("Getting nomatches for:", highcat(dt_name))

            dt$wfo_one_nomatch$dfOrigin <- dt_name
            catn("Adding an origin column with:", highcat(dt_name))
          }

          # rbind the wfo_one_nomatch data frame from each list
          combined_df <- rbind(combined_df, dt$wfo_one_nomatch)
        }
      }

      nomatch_combined <- combined_df[!duplicated(combined_df[[column]]), ]

      catn("Writing combined no-matches to:", colcat(files_dir, color = "output"))

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
