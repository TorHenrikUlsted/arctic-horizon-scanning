union_dfs <- function(df1, df2, verbose = F) {
  
  if (is.null(df1) & is.null(df2)) {
    return(data.table())
  } else if (is.null(df1) || is.null(df2)) {
    if (!is.null(df1)) {
      return(df1)
    } else if (!is.null(df2)) {
      return(df2)      
    } 
  }
  
  
  catn("Combining data tables using union.")
  
  df1_name <- deparse(substitute(df1))
  df2_name <- deparse(substitute(df2))
  
  vebcat("Merging", highcat(df1_name), "and", highcat(df2_name), veb = verbose)
  
  ## combine present lists and remove duplicates
  merged_df = dplyr::union(df1, df2)
  
  cat(sprintf("%13s | %13s | %9s \n", df1_name, df2_name, "merged_df"))
  cat(highcat(sprintf("%13d | %13d | %9d \n", nrow(df1), nrow(df2), nrow(merged_df))))
  
  catn("Duplicated species removed:", highcat(nrow(df1) + nrow(df2) - nrow(merged_df)))
  
  ## Run NA and distinct check
  if (any(is.na(merged_df)) == T) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    
    merged_df <- na.omit(merged_df)
  }
  
  if (any(duplicated(merged_df)) == T) {
    vebcat("Some merged_df species are duplicated.", color = "nonFatalError")
    
    merged_df <- unique(merged_df)
  }
  
  return(merged_df)
}

select_wfo_column <- function(filepath, col.unique, col.select = NULL, col.combine = NULL, pattern = "*.csv", verbose = FALSE) {
  
  catn("Selecting WFO column")
  
  # Check if filepath is a directory or a specific file
  if (file.info(filepath)$isdir) {
    csv_files <- list.files(path = filepath, pattern = "*.csv", full.names = TRUE)
  } else {
    csv_files <- filepath
  }
  
  # make a list of data frames based on the different CSV files and also check for any "no matches" then add those to their own data frame.
  df_list <- lapply(csv_files, function(file) {
    df <- fread(file)
    
    if (is.vector(col.unique) && length(col.unique) > 1) {
      vebcat("\nFound vector more than 1 in length, combining columns:", highcat(col.unique), veb = verbose)
      df$refinedScientificName <- apply(df[, ..col.unique, drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " "))
      df$refinedScientificName <- trimws(df$refinedScientificName)
      col.unique <- "refinedScientificName"
    }
    
    if (!is.null(col.select)) {
      vebcat("Selecting the", highcat(col.select), "column(s) and using", highcat(col.unique), "as the unique column.", veb = verbose)
      df_sel <- df %>% 
        select(all_of(c(col.select, col.unique)))
    } else {
      vebcat("Using the", col.unique, "as unique column.", veb = verbose)
      df_sel <- df %>% 
        select(all_of(col.unique))
    }
    
    df_uniq <- df_sel[!duplicated(df_sel[[col.unique]]), ]
    
    n_orig <- nrow(df_sel)
    n_uniq <- nrow(df_uniq)
   if (verbose) {
     catn("\nList:", highcat(sub("-wfo-one.csv$", "", basename(file))))
     cat(sprintf("%-10s | %s \n", "n_species", "unique(n_species)"))
     cat(highcat(sprintf("%-10d | %d \n", n_orig, n_uniq)))
     catn()
   }
    
    return(df_uniq)
  })
  
  names(df_list) <- sub("-wfo-one.csv$", "", basename(csv_files))
  names(df_list) <- gsub("-", "_", names(df_list))
  
  
  return(df_list = df_list)
}

fix_nomatches <- function(dfs, nomatch.edited, column, verbose = FALSE) {
  catn("Fixing nomatches.")
  
  # Combine no-matches
  edited_nomatch <- fread(nomatch.edited, encoding = "UTF-8")
  
  # Remove NAs and blank values in the newName df
  formatted <- edited_nomatch %>% 
    dplyr::filter(!is.na(acceptedName) & acceptedName != "") %>% 
    dplyr::select(acceptedName, listOrigin)
  
  vebprint(formatted, text = "Manually formatted data:", veb = verbose)
  
  for(i in 1:nrow(formatted)){
    
    # Get the scientificName and listOrigin from the current row
    scientificName <- formatted[i, acceptedName]
    listOrigin <- as.character(formatted[i, "listOrigin"])
    
    vebcat("Appending", highcat(as.character(scientificName)), "to", highcat(listOrigin), veb = verbose)
    
    # Check if scientificName or listOrigin is missing
    if (is.na(scientificName) | is.na(listOrigin)) {
      catn("Missing value in row ", i)
      next
    }
    
    # Check if listOrigin exists in dfs
    if (!listOrigin %in% names(dfs)) {
      vebcat("Invalid listOrigin in row ", i, ": ", listOrigin, color = "nonFatalError")
      next
    }
    
    # Create a new data table with the same columns as dfs[[listOrigin]]
    new_row <- data.table(scientificName)
    setnames(new_row, names(dfs[[listOrigin]]))
    
    # Append the scientificName to the corresponding data table in dfs
    dfs[[listOrigin]] <- rbindlist(list(dfs[[listOrigin]], new_row))
  }
  
  vebcat("the manually formatted synonym checks have been successfully added to correct data frames.", veb = verbose, color = "proSuccess")
  
  return(dfs)
}

write_filter_fun <- function(file.out, spec.in, fun) {
  create_dir_if(dirname(file.out))
  
  catn(highcat(nrow(spec.in)), "Species input.")
  
  result <- fun()
  
  catn(highcat(nrow(spec.in) - nrow(result)), "Species removed:")
  
  catn(highcat(nrow(result)), "species output")
  
  catn("Writing species to:", colcat(file.out, color = "output"))
  
  fwrite(result, file.out, bom = T)
  
  return(result)
}

# For file.out, do not use extension
get_occurrence <- function(spec, file.out, region = NULL, coord.uncertainty = 4600, download.key = NULL, download.doi = NULL, verbose = FALSE) {
  
  spec_keys <- get_sp_keys(
    sp_names = spec,
    out.dir = dirname(file.out),
    verbose = verbose
  )
  
  occ_data <- get_occ_data(
    species_w_keys = spec_keys,
    file.name = file.out,
    region = region,
    coord.uncertainty = coord.uncertainty,
    download.key = download.key,
    download.doi = download.doi,
    verbose = verbose
  )
  
  return(list(
   occ = occ_data,
   keys = spec_keys
  ))
}

get_prefix <- function(taxonRank) {
  switch(taxonRank,
         "SPECIES" = "",
         "SUBSPECIES" = "subsp.",
         "VARIETY" = "var.",
         "FORM" = "f.",
         "")
  
}

combine_columns <- function(dt, col1, col2, col3, verbose = FALSE) {
  if (is.data.table(dt)) {
    
    catn("Creating prefix column.")
    
    dt[, prefix := sapply(dt[[col3]], get_prefix)]
    
    dt[, combined := trimws(do.call(paste, c(.SD, sep = " "))), .SDcols = c(col1, "prefix", col2)]
    
  } else {
    stop("Error: Input is not a data table.")
  }
  
  return(dt)
}

remove_authorship <- function(dt, verbose = FALSE) {
  vebcat("Removing Authorship from scientificName.", veb = verbose)
  
  signs <- c("×")
  
  dt[, cleanName := sapply(scientificName, function(x) {
    components <- strsplit(x, " ")[[1]]
    result <- c(components[1], components[2])
    
    if (components[2] %in% signs) {
      result <- c(result, components[3])
    }
    
    if (components[3] %in% infraEpithet_designations) {
      result <- c(result, components[3], components[4])
    }
    
    if (components[4] %in% infraEpithet_designations || components[4] %in% signs) {
      result <- c(result, components[4], components[5])
    }
    
    if (components[5] %in% signs) {
      result <- c(result, components[5], components[6])
    }
    
    
    return(paste(result, collapse = " "))
  })]
  
  return(dt)
}

# spec.occ can be either filepath or data frame
chunk_protocol <- function(
    spec.occ,
    spec.keys,
    chunk.name = "species",
    chunk.col,
    chunk.dir,
    chunk.size = 1e6,
    cores.max = 1,
    iterations = NULL,
    verbose = FALSE
  ) {
  
  vebprint(head(spec.occ, 1), text = "spec.occ", veb = verbose)
  
  if (is.character(spec.occ)) {
    chunk_file(
      file_path = spec.occ,
      chunk.name = chunk.name,
      chunk.column = chunk.col, 
      chunk.dir = chunk.dir, 
      chunk.size = chunk.size,
      cores.max = cores.max,
      iterations = iterations,
      verbose = verbose
    )
    
  } else if ("data.frame" %in% class(spec.occ) || "data.table" %in% class(spec.occ) ) {
    chunk_loaded_df(
      df = spec.occ,
      chunk.name = chunk.name,
      chunk.column = chunk.col, 
      chunk.dir = chunk.dir, 
      verbose = verbose
    )
    
  } else {
    stop("Invalid input, either filepath or data frame/table")
  }
  
  # If chunk_loaded_df or chunk_file has a vector input, then "combined" must be used as chunk.column parameter, else the same as chunk_loaded_df. chunk.name has to be the same for all.
  if (is.vector(chunk.col)) {
    chunk.col = "cleanName"
  }

  clean_chunks(
    chunk.name = chunk.name,
    chunk.column = chunk.col, 
    chunk.dir = chunk.dir, 
    sp_w_keys = spec.keys,
    verbose = verbose
  )
}


