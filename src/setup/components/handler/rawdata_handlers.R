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
  }

  vebcat("List name:", highcat(name), veb = verbose)
  vebcat("Directory:", highcat(dir), veb = verbose)

  create_dir_if(dir)
  present_out <- paste0(dir, "/", name, "-present.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")

  if (!file.exists(present_out)) present_out <- NULL
  if (!file.exists(absent_out)) absent_out <- NULL

  if ((!is.null(present_out) || !is.null(absent_out))) {
    vebcat(highcat(name), "already wrangled, loading files..", veb = verbose)

    res <- list()

    if (!is.null(present_out)) res$present <- fread(present_out)
    if (!is.null(absent_out)) res$absent <- fread(absent_out)
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

      vebprint(wrangle_funs, verbose, text = "all Wrangle functions found:")

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
      vebcat("Error occurred when trying to wrangle all files.", color = "fatalError")
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

    check_file <- paste0(out_dir_child, "/wfo-completed.txt")

    file_name <- paste0(parent_folder, "-", child_folder, "-wfo")

    if (!file.exists(check_file)) {
      # Run synonym check on the species
      sp_synonyms <- check_syn_wfo(
        checklist = wrangled_dfs[[name]],
        cols = list(
          spec.name = column,
          Authorship = if(paste0(column, "Authorship") %in% names(wrangled_dfs[[name]])) paste0(column, "Authorship")
        ),
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

      manual_checks <- sum(
        nrow(sp_synonyms$mismatch),
        nrow(sp_checked$nomatch),
        nrow(sp_checked$na),
        na.rm = TRUE
      )

      lost_diff <- (
        nrow(wrangled_dfs[[name]]) - (
          nrow(sp_checked$clean) + manual_checks + nrow(sp_checked$duplicate)
        )
      )
      
      if (lost_diff < 0) lost_diff <- 0

      md_dt <- data.table(
        wrangle = nrow(wrangled_dfs[[name]]),
        match = nrow(sp_synonyms$clean) + nrow(sp_synonyms$mismatch),
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

      sp_checked$raw <- NULL
      
      sp_checked$mismatch <- rbind(
        unique(sp_checked$mismatch, by = column),
        unique(sp_synonyms$mismatch, by = column),
        fill = TRUE
      )
      
      # Choose approach
      if (config$simulation$approach == "precautionary") {
        wfo_data <- WFO.minimal(WFO_file)
        
        sp_checked$clean <- filter_approach(
          dt = sp_checked$clean,
          WFO.data = wfo_data,
          out.file = paste0(out_dir_child, "/wfo-one-approach.csv"),
          verbose
        )

        sp_checked$mismatch <- filter_approach(
          dt = sp_checked$mismatch,
          WFO.data = wfo_data,
          out.file = paste0(out_dir_child, "/wfo-one-approach-mismatch.csv"),
          verbose
        )
      }

      create_file_if(check_file)

      return(setNames(list(sp_checked), name))
    } else {
      vebcat(highcat(name), "already synonym checked.", veb = verbose)
    }
  })

  synonym_lists <- unlist(synonym_lists, recursive = FALSE)

  return(synonym_lists)
}

setup_raw_data <- function(column, cores.max = 1, verbose = FALSE, counter = 1) {
  vebcat("Setting up raw data.", color = "funInit")

  files_dir <- "./outputs/setup/wrangle"

  create_dir_if(files_dir)
  create_dir_if(paste0(files_dir, "/test"))

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
      vebcat("Combining no-matches.", color = "proInit")
      
      test_files_dir <- paste0(files_dir, "/test")
      
      combined_dt_test <- data.table()
      combined_dt <- data.table()
      
      # Loop over each list in checked_dfs
      for (dt_name in names(checked_dts)) {
        wfo_one <- checked_dts[[dt_name]]
        
        if (!is.null(wfo_one$clean)) wfo_one$clean <- NULL
        if (!is.null(wfo_one$duplicate)) wfo_one$duplicate <- NULL
        
        vebprint(wfo_one, verbose, paste0("Checking data table ", dt_name, " :"))
        process_status <- function(status) {
          if (!is.null(wfo_one[[status]]) && nrow(wfo_one[[status]]) > 0) {
            catn("Appending", highcat(nrow(wfo_one[[status]])), status, "for:", highcat(dt_name))
            wfo_one[[status]][, `:=`(listOrigin = dt_name, status = status)]
            return(wfo_one[[status]])
          }
          return(NULL)
        }
        
        status <- c("mismatch", "nomatch", "na")
        processed <- lapply(status, process_status)
        total_rows <- sum(sapply(processed, function(x) if (!is.null(x)) nrow(x) else 0))
        
        if (total_rows > 0) {
          is_test <- grepl("test", dt_name)
          current_dt <- if (is_test) combined_dt_test else combined_dt
          
          current_dt <- rbindlist(c(list(current_dt), processed), fill = TRUE)
          
          if (is_test) {
            combined_dt_test <- current_dt
          } else {
            combined_dt <- current_dt
          }
        }
        
        mdwrite(
          config$files$post_seq_md,
          text = paste0("3;Combined **", dt_name, "** for manual handling for: **", total_rows, "**")
        )
      }
      
      process_manual_check <- function(dt, column, out.dir, is.test) {
        if (nrow(dt) > 0) {
          file_mismatch <- paste0(out.dir, "/combined-wfo-mismatches.csv")
          file_manual <- paste0(out.dir, "/manual-check-file.csv")
          # After the loop, process the combined data
          dt <- unique(dt, by = column)
          
          catn("There were", highcat(nrow(dt)), 
               ifelse(is.test, "test", ""), 
               "species that need manual handling"
              )
          
          catn("Writing manual edit combined to:", colcat(file_mismatch, color = "output"))
          fwrite(dt, file_mismatch, bom = TRUE)
          
          man_dt <- data.table(
            verbatimName = dt$verbatimName,
            interimName = dt[[column]],
            interimNameAuthorship = if(paste0(column, "Authorship") %in% names(dt)) dt[[paste0(column, "Authorship")]],
            wfoName = dt$scientificName.ORIG,
            wfoSpecies = dt$scientificName,
            Old.status = dt$Old.status,
            Old.name = dt$Old.name,
            name.clean = dt$name.clean,
            mismatch.old = dt$mismatch.old,
            mismatch.scientific = dt$mismatch.scientific,
            status = dt$status,
            acceptedName = NA_character_,
            acceptedNameAuthorship = NA_character_,
            source = NA_character_,
            comment = NA_character_,
            listOrigin = dt$listOrigin
          )
          
          catn("Use this file to manually edit:", colcat(file_manual, color = "output"))
          # Add help link
          fwrite(man_dt, file_manual, bom = TRUE)
          
          vebcat("no-matches combined successfully.", color = "proSuccess")
        } else {
          catn("There were", highcat(0), ifelse(is.test, "test", ""), "species in need of manual handling \n")
        }
      }
      
      # Process test and non-test data
      if (nrow(combined_dt_test) > 0) {
        process_manual_check(combined_dt_test, column, test_files_dir, TRUE)
      }
      process_manual_check(combined_dt, column, files_dir, FALSE)
      
    }
  }

  vebcat("raw data setup completed successfully.", color = "funSuccess")
}
