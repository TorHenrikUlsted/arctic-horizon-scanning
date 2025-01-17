#------------------------#
####      System      ####
#------------------------#

source_all <- function(dir) {
  # Get a list of all .R files in the directory and its subdirectories
  r_files <- list.files(path = dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
  
  # Source each file
  lapply(r_files, source)
  
  cat(length(r_files), "scripts sourced from", dir, "and its subdirectories.\n")
}

system_calc_rows <- function(file.path) {
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    # Use PowerShell to count lines
    total_rows <- as.numeric(system2("powershell", 
                                     args = c("-command", 
                                              sprintf("(Get-Content '%s' | Measure-Object -Line).Lines", 
                                                      file.path)), 
                                     stdout = TRUE))
  } else {
    # Unix-based systems
    total_rows <- as.numeric(system(paste("awk 'END {print NR}' ", file.path), intern = TRUE))
  }
  # -1 for header
  return(total_rows - 1)
}

system_calc_unique <- function(file.path, column, sep = "\t") {
  os <- Sys.info()["sysname"]
  
  if (is.character(column)) {
    tmp <- fread(file.path, nrows = 0)
    column_number <- which(names(tmp) == column)
    if (length(column_number) == 0) stop("Column name not found")
  } else if (is.numeric(column)) {
    column_number <- as.integer(column)
  } else {
    stop("Column must be either a name (character) or a number")
  }
  
  if (os == "Windows") {
    # Windows approach using R
    unique_values <- new.env(hash = TRUE)
    con <- file(file.path, "r")
    header <- readLines(con, n = 1)  # Skip header
    while (length(line <- readLines(con, n = 1)) > 0) {
      fields <- strsplit(line, sep, fixed = TRUE)[[1]]
      if (column_number <= length(fields)) {
        unique_values[[fields[column_number]]] <- TRUE
      }
    }
    close(con)
    unique_count <- length(unique_values)
  } else {
    # Unix-based systems
    awk_sep <- gsub("\\", "\\\\", sep, fixed = TRUE)
    cmd <- sprintf("awk -F'%s' '{print $%d}' '%s' | sort -u | wc -l", 
                   awk_sep, column_number, file.path)
    unique_count <- as.numeric(system(cmd, intern = TRUE))
  }
  # -1 for header
  return(unique_count - 1)
}

system_calc_uniq_and_rows <- function(file.path, column = NULL, sep = "\t") {
  if (is.null(column)) {
    stop("Column must be specified")
  }
  
  if (is.character(column)) {
    tmp <- fread(file.path, nrows = 0)
    column_number <- which(names(tmp) == column)
    if (length(column_number) == 0) stop("Column name not found")
  } else if (is.numeric(column)) {
    column_number <- as.integer(column)
  } else {
    stop("Column must be either a name (character) or a number")
  }
  
  os <- Sys.info()["sysname"]
  
  if (os == "Windows") {
    # Windows approach using PowerShell
    ps_script <- sprintf(
      "Get-Content '%s' | Select-Object -Skip 1 | ForEach-Object {
         $_.Split('%s')[%d]
       } | Group-Object | Measure-Object Count,Name",
      file.path, sep, column_number - 1
    )
    result <- shell(sprintf("powershell -Command \"%s\"", ps_script), intern = TRUE)
    
    # Parse PowerShell output (-1 for header)
    total_rows <- as.numeric(strsplit(result[1], "\\s+")[[1]][2]) - 1
    unique_count <- as.numeric(strsplit(result[2], "\\s+")[[1]][2]) - 1
  } else {
    # Unix-based systems
    awk_sep <- gsub("\\", "\\\\", sep, fixed = TRUE)
    cmd <- sprintf("awk -F'%s' '{print $%d}' '%s' | awk '{count++; unique[$0]++} END {print count, length(unique)}'",
                   awk_sep, column_number, file.path)
    result <- system(cmd, intern = TRUE)
    values <- as.numeric(strsplit(result, " ")[[1]])
    total_rows <- values[1] - 1  # -1 for header
    unique_count <- values[2] - 1  # -1 for header
  }
  
  return(list(
    total_rows = total_rows,
    unique_count = unique_count
  ))
}

create_derived_dataset <- function(data.name, verbose = FALSE) {
  if (data.name$spec.unknown == "" || is.na(data.name$spec.unknown)) {
    stop("Need to include a name for at least spec.unknown")
  }
  
  run <- 1
  
  filter_dir <- "./outputs/filter/"
  chunk_dir <- "/chunk/species"
  post_process_dir <- "./outputs/post-process"
  
  repeat {
    if (config$simulation$validation && run == 2 && data.name$spec.known != "" || is.na(data.name$spec.known)) {
      derived_data <- paste0("derived-data-", data.name$spec.known)
      occ_dir <- paste0(filter_dir, data.name$spec.known, chunk_dir)
      sp_occ_out <- file.path(post_process_dir, derived_data, "datasetKey-count.csv")
      derived_data_zip_out <- file.path(post_process_dir, derived_data, "derived-dataset.zip")
    } else {
      derived_data <- paste0("derived-data-", data.name$spec.unknown)
      occ_dir <- paste0(filter_dir, data.name$spec.unknown, chunk_dir)
      sp_occ_out <- file.path(post_process_dir, derived_data, "datasetKey-count.csv")
      derived_data_zip_out <- file.path(post_process_dir, derived_data, "derived-dataset.zip")
    }
    
    create_dir_if(dirname(sp_occ_out))
    
    if (file.exists(sp_occ_out)) {
      catn("DatasetKey with occurrence count already exists.")
      
      sp_occ_dt <- fread(sp_occ_out)
      
      catn(highcat(nrow(sp_occ_dt)), "datasetKeys added to csv file with corresponding total count.")
      
      catn(highcat(sum(sp_occ_dt$count)), "Total occurrences.")
      
      rm(sp_occ_dt)
      invisible(gc())
      
      if (!config$simulation$validation || run >= 2) break
      run <- run + 1
    } else {
      vebcat("Collecting datasetKeys and corresponding occurrence counts", color = "funInit")
      
      sp_occ_files <- list.files(occ_dir, full.names = TRUE)
      
      vebprint(head(sp_occ_files, 5), verbose, "occurence files:")
      
      sp_occ_dt <- data.table(datasetKey = character(), count = integer())
      
      catn("Reading data.")
      for (i in 1:length(sp_occ_files)) {
        sp_occ <- sp_occ_files[i]
        
        vebprint(sp_occ, verbose, "sp_occ:")
        
        cat("\rCounting datasetKey occurrences for", i, "/", length(sp_occ_files))
        
        sp_dt <- fread(sp_occ, select = "datasetKey")
        
        vebprint(sp_dt, verbose, "sp_dt:")
        
        sp_count <- sp_dt[, .(count = .N), by = datasetKey]
        
        vebprint(sp_count, verbose, "sp_count:")
        
        sp_occ_dt <- rbind(sp_occ_dt, sp_count)
      }
      catn()
      
      sp_occ_dt <- sp_occ_dt[, .(count = sum(count)), by = datasetKey]
      
      if (any(duplicated(sp_occ_dt$datasetKey))) {
        vebcat("Some datasetKeys are duplicated.", color = "nonFatalError")
      } else {
        vebcat("All datasetKeys are unique.", color = "proSuccess")
      }
      
      catn(highcat(nrow(sp_occ_dt)), "datasetKeys added to csv file with corresponding total count.")
      
      catn(highcat(sum(sp_occ_dt$count)), "Total occurrences.")
      
      vebprint(sp_occ_dt, text = "Sample output:")
      
      catn("Writing out to file:", colcat(sp_occ_out, color = "output"))
      fwrite(sp_occ_dt, sp_occ_out)
      
      vebcat("Successfully collected datasetKeys and counts", color = "funSuccess")
    }
    
    if (file.exists(derived_data_zip_out)) {
      catn("Derived dataset already exists, skipping process.")
      files <- list.files(occ_dir, full.names = FALSE)
      vebcat("Found", highcat(length(files)), "species in the directory.")
      rm(files)
      invisible(gc())
    } else {
      files <- list.files(occ_dir, full.names = TRUE)
      
      catn("Zipping", highcat(length(files)), "dervied species files.")
      
      zip(
        zipfile = derived_data_zip_out, 
        files = files,
        flags = "-j" # remove directories and only keep files
      )
    }
    
    if (!config$simulation$validation || run >= 2) break
    run <- run + 1
  }
}

progressive_dirname <- function(path, begin = 1, end = NULL) {
  # Split the path into directories
  dirs <- strsplit(path, "/")[[1]]
  
  # If end is NULL, include all directories from begin to the end of the path
  if (is.null(end)) {
    end <- length(dirs)
  }
  
  # Subset the directories based on begin and end
  dirs_subset <- dirs[begin:end]
  
  # Combine the subset of directories back into a path
  res <- paste(dirs_subset, collapse = "/")
  
  return(res)
}

remove_parent_paths <- function(paths, verbose = FALSE) {
  # Sort paths by length to optimize comparisons
  # (longer paths might be nested under shorter ones)
  paths <- sort(paths, decreasing = FALSE)
  
  vebprint(paths, verbose, "Input paths (sorted):")
  
  # For each path, check if it's a prefix of any other path
  is_parent <- sapply(paths, function(potential_parent) {
    # Add trailing separator to ensure we match complete path components
    parent_with_sep <- paste0(potential_parent, .Platform$file.sep)
    # Check if this path is a prefix of any other path
    any(startsWith(paths, parent_with_sep))
  })
  
  # Keep only paths that aren't parents of other paths
  result <- paths[!is_parent]
  
  if(verbose) {
    cat("\nRemoved parent paths:\n")
    print(paths[is_parent])
    cat("\nFinal paths:\n")
    print(result)
  }
  
  return(result)
}

truncate_vector <- function(x, length.max = 20, verbose = FALSE) {
  if (is.null(x)) return(invisible())
  if (length(x) <= length.max) return(invisible(x))
  
    if (is.character(x) && any(!is.na(file.info(x)$isdir))) {
      if (!all(file.info(x)$is.dir)) { # some are files
        x <- head(x, length.max)
        vebcat("Truncated with head for exceeding length of", length.max, veb = verbose)
      } else {
        x <- unique(dirname(x))
        vebcat("Truncated by dirname for exceeding length of", length.max, veb = verbose)
      }
      
    } else {
      stop("Wrong class: ", class(x))
    }
  
  return(x)
}

get_files <- function(input.dir, exclude.dirs = NULL, exclude.files = NULL, include.files = NULL, step = 0, verbose = FALSE) {
  if (step > 6) {
    vebcat("Step cannot be higher than 6", color = "fatalError")
    stop("Change step to a different integer value.")
  }
  max_length <- 20
  # Get all directories
  d_orig <- list.dirs(input.dir, recursive = TRUE)
  if (length(d_orig) > 1) {
    d_orig <- d_orig[-1]
  } else {
    stop("Could not find any directories.")
  }
  d_orig_u <- truncate_vector(d_orig, max_length, step == 1)
    vebprint(
      remove_parent_paths(d_orig_u), veb = (step == 1 || verbose), 
      paste0(
        highcat("Step 1"), " | ",
        highcat("initial directories [", length(d_orig), "] :")
      )
    )
  if (step == 1) return(invisible())
  
  # Remove exclude.dirs
  if (!is.null(exclude.dirs)) {
    exclude <- sapply(d_orig, function(dir) any(sapply(exclude.dirs, function(ex_d) grepl(ex_d, dir))))
    d <- d_orig[!exclude]
  }
  d_u <- truncate_vector(d, max_length, step == 2)
  vebprint(
    d_u, veb = (step == 2 || verbose), 
    paste0(
      highcat("Step 2"), " | ", 
      highcat("After removing exclude.dirs [", length(d),"] :")
    )
  )
  if (step == 2) return(invisible())
  
  # List all files and find the include.files
  if (!is.null(include.files)) {
    f_inc <- list.files(d_orig, recursive = FALSE, full.names = TRUE)
    include <- sapply(f_inc, function(file) any(sapply(include.files, function(inc_f) grepl(inc_f, file, fixed = TRUE))))
    f_inc <- f_inc[include]
    f_inc <- remove_parent_paths(f_inc)
  } else {
    f_inc <- NULL
  }
  f_inc_u <- truncate_vector(f_inc, max_length, step == 3)
  vebprint(
    f_inc_u, veb = (step == 3 || verbose), 
    paste0(
      highcat("Step 3"), " | ", 
      highcat("Files included in include.files [", length(f_inc), "] :")
    )
  )
  if (step == 3) return(invisible())
  
  # Get the files in the included directories
  f <- list.files(d, recursive = FALSE, full.names = TRUE)
  f <- f[!file.info(f)$isdir]
  f_u <- truncate_vector(f, max_length, step == 4)
  vebprint(
    f_u, veb = (step == 4 || verbose), 
    paste0(
      highcat("Step 4"), " | ", 
      highcat("Files in the included directories [", length(f), "] :")
    )
  )
  if (step == 4) return(invisible())
  
  # Remove the exclude.files
  if (!is.null(exclude.files)) {
    exclude <- sapply(f, function(x) any(sapply(exclude.files, function(ex_f) grepl(ex_f, x))))
    f <- f[!exclude]
    f <- remove_parent_paths(f)
    f <- c(f, f_inc)
  }
  f_u <- truncate_vector(f, max_length, step == 5)
  vebprint(
    f_u, veb = (step == 5 || verbose), 
    paste0(
      highcat("Step 5"), " | ", 
      highcat("After removing exclude.files [", length(f), "] :")
    )
  )
  if (step == 5) return(invisible())
  
  return(f)
}

get_repository_files <- function(which.sequence = "all", step = 0, subset = NULL, verbose = FALSE) {
  vebcat("Collecting repository files", color = "funInit")
  
  dirs <- list(
    setup = list(
      dir = "./outputs/setup",
      exclude_dirs = c("wfo-match-nodes", "test", "logs", "system", "locks", "projections", "region", "old"),
      exclude_files = c("stop-file.txt", ".zip", "wfo-match-nodes", ".rds", ".tif", "stats.csv", "wfo-completed"),
      include_files = c("coordinateUncertainty"),
      step = step
    ),
    filter = list(
      dir = "./outputs/filter",
      exclude_dirs = c("chunk", "test", "old"),
      exclude_files = c("occ.csv", ".zip", "logs", "completed.txt"),
      step = step
    ),
    hypervolume = list(
      dir = "./outputs/hypervolume",
      exclude_dirs = c("test", "prep", "locks", "old"),
      exclude_files = c("init-log.txt", "node-iterations.txt"),
      step = step
    ),
    visualize = list(
      dir = "./outputs/visualize",
      exclude_dirs = c("stack", "test"),
      exclude_files = c("warning", "error"),
      step = step
    ),
    post_process = list(
      dir = "./outputs/post-process",
      exclude_dirs = c("old"),
      step = step
    ),
    utils = list(
      dir = "./outputs/utils",
      exclude_dirs = "logs",
      step = step
    )
  )
  
  repo_files <- c()
  
  for (sequence in names(dirs)) {
    if (which.sequence == "all" || which.sequence == sequence) {
      catn("Collecting", sequence, "files")
      if (verbose) print_function_args()
      
      tryCatch({
        files <- get_files(
          input.dir = dirs[[sequence]]$dir,
          exclude.dirs = dirs[[sequence]]$exclude_dirs,
          exclude.files = dirs[[sequence]]$exclude_files,
          include.files = dirs[[sequence]]$include_files,
          step = dirs[[sequence]]$step,
          verbose = verbose
        )
      }, error = function(e) {
        stop(e)
      })
      
      if (!is.null(subset)) {
        subset_dir <- paste0(dirs[[sequence]]$dir, "/", subset)
        
        indices <- grep(paste0("^", subset_dir), files)
        
        files <- files[indices]
      }
      
      repo_files <- c(repo_files, files)
    }
  }
  
  
  vebcat("Repository files collected succesfully", color = "funSuccess")
  
  if (step == 0) return(repo_files) else return(invisible(repo_files))
}

pack_repository <- function(filename = "Horizon-Scanning-Repository", which.sequence = "all") {
  vebcat("Packing repository", color = "funInit")
  
  out_file <- paste0("./outputs/", filename, ".zip")
  
  if (file.exists(out_file)) {
    catn("found file:", colcat(out_file, color = "output"))
    vebcat("Repository already zipped, remove it or rename the file.", color = "fatalError")
    return(invisible())
  }
  
  repo_files <- get_repository_files(which.sequence = which.sequence)
  
  catn("Compressing files with .zip")
  zip(zipfile = out_file, files = repo_files)
  
  catn("Zip file saved in", colcat(out_file, color = "output"))
  
  vebcat("Repository packed successfully", color = "funSuccess")
}

create_clickable_link <- function(url, text = NULL, file = TRUE, params = NULL) {
  if (is.null(text)) text <- as.character(url)
  url <- normalizePath(url, winslash = "/")
  
  if (file) {
    link <- cli::cli_text("{.href [{text}](file://{url})}")
  } else if (!file && cli::ansi_has_hyperlink_support()) {
    link <- cli::style_hyperlink(text, url, params = NULL)
  }
  
  return(link)
}

find_term <- function(term, dir = ".", file.pattern = "\\.R$", file.exclude = NULL, verbose = FALSE) {
  term <- to_char(term, verbose = verbose)
  files <- list.files(dir, pattern = file.pattern, full.names = TRUE, recursive = TRUE)
  
  if (!is.null(file.exclude)) {
    files <- files[!sapply(files, function(file) {
      any(sapply(file.exclude, function(exclude) grepl(exclude, file, fixed = TRUE)))
    })]
  }
  
  results <- lapply(files, function(file) {
    tryCatch(
      {
        lines <- readLines(file, warn = FALSE)
        # Use word boundaries and lookahead/lookbehind for exact matches
        pattern <- paste0("(?<!(\\w|\\$))", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", term), "(?!(\\w|\\$))")
        matches <- which(sapply(lines, function(line) grepl(pattern, line, perl = TRUE)))
        
        if (length(matches) > 0) {
          if (verbose) {
            #cat("File:", file, "\n")
            cat("Matches found on lines:", paste(matches, collapse = ", "), "\n")
            cat("Matching lines:\n")
            for (m in matches) {
              cat("  Line", m, ":", lines[m], "\n")
            }
            cat("\n")
          }
          
          link_text <- create_clickable_link(file, basename(file))
          cat(link_text)
          
          data.table(
            file = file,
            lineNumber = matches,
            code = trimws(lines[matches]),
            stringsAsFactors = FALSE
          )
        } else {
          NULL
        }
      },
      error = function(e) {
        if (verbose) warning("Error processing file: ", file, "\nError: ", e$message)
        NULL
      }
    )
  })
  return(rbindlist(results[!sapply(results, is.null)]))
}

find_term_pattern <- function(term, line.pattern = NULL, file.exclude = NULL, dir = ".", file.pattern = "\\.R$") {
  term <- to_char(term)
  line.pattern <- to_char(line.pattern)
  
  find_res <- find_term(term, dir, file.pattern, file.exclude)
  
  if (is.null(find_res)) {
    message("No matches found.")
    return(NULL)
  }
  
  if (!is.null(line.pattern)) {
    # General case for function calls
    if (grepl("\\(\\)$", line.pattern)) {
      base_pattern <- sub("\\(\\)$", "", line.pattern)
      line.pattern <- paste0(base_pattern, "\\([^\\)]*\\)")
    } else if (grepl("\\[\\]$|\\{\\}$", line.pattern)) {
      base_pattern <- sub("\\[\\]$|\\{\\}$", "", line.pattern)
      bracket_type <- substring(line.pattern, nchar(line.pattern) - 1, nchar(line.pattern))
      
      escape_chars <- list("[]" = "\\[\\]", "{}" = "\\{\\}")
      escaped_bracket <- escape_chars[[bracket_type]]
      
      line.pattern <- paste0(base_pattern, escaped_bracket[1], "[^", escaped_bracket[2], "]*", escaped_bracket[2])
    }
    
    
    pattern_matches <- grepl(line.pattern, find_res$code)
    filtered_res <- find_res[pattern_matches, ]
    
    removed_count <- nrow(find_res) - nrow(filtered_res)
    cat("Number of results without the pattern:", removed_count, "\n")
    
    if (nrow(filtered_res) == 0) {
      message("No matches found after applying the line pattern.")
      return(NULL)
    }
    
    return(filtered_res)
  } else {
    return(find_res)
  }
}

remove_outer_pattern <- function(text, pattern, replacement = "", show.diff = TRUE, verbose = FALSE) {
  # Split the pattern into start and end
  split_pattern <- function(p) {
    if (p == "") {
      return(list(start = "", end = ""))
    }
    parts <- strsplit(p, "")[[1]]
    open_bracket <- which(parts %in% c("(", "[", "{"))
    close_bracket <- which(parts %in% c(")", "]", "}"))
    if (length(open_bracket) == 0 || length(close_bracket) == 0) {
      stop("Invalid pattern: must contain opening and closing brackets")
    }
    list(
      start = paste(parts[1:open_bracket], collapse = ""),
      end = paste(parts[close_bracket:length(parts)], collapse = "")
    )
  }
  
  pattern_parts <- split_pattern(pattern)
  replacement_parts <- split_pattern(replacement)
  
  vebprint(pattern_parts, verbose, "Pattern parts:")
  vebprint(replacement_parts, verbose, "Replacement parts:")
  
  # Escape special characters in the patterns
  pattern_start_escaped <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", pattern_parts$start)
  pattern_end_escaped <- gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", pattern_parts$end)
  
  # Create the regex pattern
  full_pattern <- paste0(pattern_start_escaped, "(.*?)", pattern_end_escaped)
  
  # Replace the pattern with its contents and the specified replacement
  result <- gsub(full_pattern, paste0(replacement_parts$start, "\\1", replacement_parts$end), text)
  
  if (show.diff) {
    # Highlight the changes
    highlighted_original <- gsub(
      full_pattern,
      paste0(red(pattern_parts$start), "\\1", red(pattern_parts$end)),
      text
    )
    highlighted_result <- gsub(
      full_pattern,
      paste0(green(replacement_parts$start), "\\1", green(replacement_parts$end)),
      text
    )
    
    vebcat("Original:\n", highlighted_original, veb = verbose)
    vebcat("Result:\n", highlighted_result, veb = verbose)
    
    return(list(result = result, high_orig = highlighted_original, high_res = highlighted_result))
  } else {
    return(list(result))
  }
}

replace_term_name <- function(name.old, name.new, dir = ".", file.pattern = "\\.R$", file.exclude = NULL, verbose = FALSE) {
  name.old <- to_char(name.old, string = "Old name after check:", verbose = verbose)
  name.new <- to_char(name.new, string = "New name after check:", verbose = verbose)
  res <- find_term(name.old, dir, file.pattern, file.exclude = file.exclude, verbose = verbose)
  
  if (is.null(res)) {
    message("No occurrences of '", name.old, "' found. No changes made.")
    return(invisible(NULL))
  }
  unique_files <- unique(res$file)
  
  for (file in unique_files) {
    tryCatch(
      {
        original_lines <- readLines(file, warn = FALSE)
        pattern <- gsub("\\$", "\\\\$", name.old)
        pattern <- paste0("(^|[^[:alnum:]_])(", pattern, ")([^[:alnum:]_]|$)")
        
        new_lines <- gsub(pattern, paste0("\\1", name.new, "\\3"), original_lines)
        
        # Identify lines where the replacement occurred
        replaced_lines <- which(original_lines != new_lines)
        
        vebprint(replaced_lines, verbose, "Replaced lines:")
        
        if (length(replaced_lines) > 0) {
          # Only style if changes were made
          temp_file <- tempfile(fileext = ".R")
          writeLines(new_lines, temp_file)
          tryCatch(
            {
              styler::style_file(temp_file)
              styled_lines <- readLines(temp_file, warn = FALSE)
              writeLines(styled_lines, file)
            },
            error = function(e) {
              warning("Error during styling: ", e$message, ". Writing unstyled changes.")
              writeLines(new_lines, file)
            }
          )
          catn("Updated", file, "at line\n-", paste(replaced_lines, collapse = "\n- "), "\n")
        }
        
        # Clean up temporary file
        if (exists("temp_file")) file.remove(temp_file)
      },
      error = function(e) {
        warning("Error processing file: ", file, "\nError: ", e$message)
      }
    )
  }
  vebcat("Finished updating", paste0("'", highcat(name.old), "'"), "to", paste0("'", highcat(name.new), "'"), "in all files.", color = "proSuccess")
  if (verbose) {
    catn("Output:")
    find_term(name.new, dir, file.pattern, file.exclude)
  }
}

replace_term_pattern <- function(term, line.pattern = NULL, file.exclude = NULL, replacement = "", dir = ".", file.pattern = "\\.R$", replace = FALSE, verbose = FALSE) {
  term <- to_char(term)
  line.pattern <- to_char(line.pattern)
  replacement <- to_char(replacement)
  
  matches <- find_term_pattern(term, line.pattern, dir, file.pattern, file.exclude = file.exclude, verbose = verbose)
  
  if (is.null(matches)) {
    message("No matches found to replace.")
    invisible(return(NULL))
  }
  
  replace_in_file <- function(file, line_numbers, old_code, new_code, verbose = FALSE) {
    lines <- readLines(file, warn = FALSE)
    for (i in seq_along(line_numbers)) {
      lines[line_numbers[i]] <- new_code[i]
    }
    if (replace) {
      writeLines(lines, file)
    }
    return(data.table(file = file, lineNumber = line_numbers, oldCode = old_code, newCode = new_code))
  }
  
  # Group matches by file
  matches_by_file <- split(matches, matches$file)
  
  # Apply replacements
  results <- lapply(names(matches_by_file), function(file) {
    file_matches <- matches_by_file[[file]]
    old_code <- file_matches$code
    new_res <- remove_outer_pattern(text = old_code, pattern = line.pattern, replacement = replacement, verbose = verbose)
    high_orig <- new_res$high_orig
    high_res <- new_res$high_res
    new_code <- new_res$result
    replaced <- replace_in_file(file, file_matches$lineNumber, old_code, new_code, verbose = verbose)
    
    return(list(replaced = replaced, high_orig = high_orig, high_res = high_res))
  })
  
  replaced_data <- do.call(rbind, lapply(results, function(x) x$replaced))
  high_orig <- do.call(rbind, lapply(results, function(x) x$high_orig))
  high_res <- do.call(rbind, lapply(results, function(x) x$high_res))
  
  if (!replace) {
    catn(
      "\n",
      highcat("######################################"), "\n",
      highcat("#              SAFEMODE              #"), "\n",
      highcat("######################################"), "\n"
    )
    cat("These are the changes that would be made:\n")
    for (i in 1:nrow(replaced_data)) {
      catn(blue$bold(paste0(highcat("File: "), replaced_data$file[i], highcat(" Line: "), replaced_data$lineNumber[i])))
      catn(highcat("Old: "), high_orig[i])
      catn(highcat("New: "), high_res[i], "\n")
    }
    
    catn(highcat("\nTo apply these changes, run again with replace = TRUE"))
  } else {
    catn("Changes applied:")
    for (i in 1:nrow(replaced_data)) {
      catn(blue$bold(paste0(highcat("File: "), replaced_data$file[i], highcat(" Line: "), replaced_data$lineNumber[i])))
      catn(highcat("Old: "), high_orig[i])
      catn(highcat("New: "), high_res[i], "\n")
    }
  }
  
  invisible(replaced_data)
}

count_project_lines <- function(directory_path) {
  # Validate input directory path
  if (!dir.exists(directory_path)) {
    cli::cli_alert_error("Directory does not exist: {.path {directory_path}}")
    stop()
  }
  
  cli::cli_alert_info("Analyzing R files in {.path {directory_path}}")
  
  # List all R files recursively
  r_files <- list.files(
    path = directory_path,
    pattern = "\\.R$|\\.r$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(r_files) == 0) {
    cli::cli_alert_warning("No R files found in directory")
    return(0)
  }
  
  cli::cli_alert_success("Found {length(r_files)} R {cli::qty(length(r_files))} file{?s}")
  
  # Initialize counter and progress bar
  total_lines <- 0
  file_counts <- data.frame(
    filename = character(),
    lines = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Create progress bar
  cli::cli_progress_bar(
    total = length(r_files),
    format = "Processing files {cli::pb_bar} {cli::pb_percent}",
    clear = FALSE
  )
  
  # Process each file
  for (i in seq_along(r_files)) {
    file <- r_files[i]
    
    # Update progress bar
    cli::cli_progress_update()
    
    tryCatch({
      # Read file content
      lines <- readLines(file, warn = FALSE)
      
      # Remove empty lines and comments
      lines <- lines[nzchar(trimws(lines))]
      lines <- lines[!grepl("^\\s*#", lines)]
      
      # Update counters
      line_count <- length(lines)
      total_lines <- total_lines + line_count
      
      # Store individual file statistics
      file_counts <- rbind(file_counts, 
                           data.frame(
                             filename = basename(file),
                             lines = line_count
                           ))
    }, error = function(e) {
      cli::cli_alert_danger("Error processing {.file {basename(file)}}: {e$message}")
    })
  }
  
  # Clear progress bar
  cli::cli_progress_done()
  
  # Create return object with detailed information
  result <- list(
    total_lines = total_lines,
    file_count = length(r_files),
    file_details = file_counts,
    directory = directory_path
  )
  
  class(result) <- "code_count"
  
  return(result)
}

print.code_count <- function(x, ...) {
  cli::cli_h1("Code Analysis Results")
  
  cli::cli_alert_info("Directory: {.path {x$directory}}")
  cli::cli_alert_info("Total R files: {x$file_count}")
  cli::cli_alert_success("Total lines of code: {x$total_lines}")
  
  if (x$file_count > 0) {
    cli::cli_h2("Files breakdown")
    
    # Sort files by line count
    x$file_details <- x$file_details[order(-x$file_details$lines), ]
    
    # Print formatted table using cli_text
    for (i in seq_len(nrow(x$file_details))) {
      cli::cli_text("{.file {x$file_details$filename[i]}} - {x$file_details$lines[i]} lines")
    }
  }
}
