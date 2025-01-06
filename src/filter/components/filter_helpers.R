union_dts <- function(dt1, dt2, na.rm = TRUE, verbose = FALSE) {
  if (!is.data.table(dt1) || !is.data.table(dt2)) {
    stop("Both dt1 and dt2 must be data.tables")
  }
  catn("Combining data tables using union.")

  dt1_name <- deparse(substitute(dt1))
  dt2_name <- deparse(substitute(dt2))

  vebcat("Merging", highcat(dt1_name), "and", highcat(dt2_name), veb = verbose)

  ## combine present lists and remove duplicates
  merged_dt <- unique(rbind(dt1, dt2))
  
  ## Calculate duplicates using row counts
  duplicate_count <- (nrow(dt1) + nrow(dt2)) - nrow(merged_dt)
  
  ## Run NA and distinct check
  na_count <- sum(is.na(merged_dt))
  
  if (any(is.na(merged_dt)) == T) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    if (na.rm) merged_dt <- na.omit(merged_dt)
  }

  md_dt <- data.table(
    dt1 = nrow(dt1),
    dt2 = nrow(dt2),
    merged_dt = nrow(merged_dt),
    duplicates = duplicate_count,
    nas = na_count,
    added = nrow(merged_dt) - nrow(dt1)
  )
  setnames(md_dt, c("dt1", "dt2"), c(dt1_name, dt2_name))

  md_dt <- kable(md_dt)

  vebprint(md_dt, text = "Union summary:")

  mdwrite(
    config$files$post_seq_md,
    text = paste0("Combining data from **", dt1_name, "** with **", dt2_name, "**"),
    data = md_dt
  )
  
  return(merged_dt)
}

anti_union <- function(remove.from, remove.with, remove.by = NULL, na.rm = TRUE, verbose = FALSE) {
  catn("Anti joining data tables")
  
  dt1_name <- deparse(substitute(remove.from))
  dt2_name <- deparse(substitute(remove.with))
  
  vebcat("Merging", highcat(dt1_name), "and", highcat(dt2_name), veb = verbose)
  
  ## combine present lists and remove duplicates
  if (!is.null(remove.by)) {
    merged_dt <- remove.from[!remove.with, on = remove.by]
  } else {
    merged_dt <- remove.from[!remove.with]
  }
  
  ## Run NA and duplicate check
  na_count <- sum(is.na(merged_dt))
  
  if (any(is.na(merged_dt)) == TRUE) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    if (na.rm) merged_dt <- na.omit(merged_dt)
  }
  
  duplicate_count <- sum(duplicated(merged_dt))
  
  if (any(duplicated(merged_dt)) == TRUE) {
    merged_dt <- unique(merged_dt)
  }
  
  
  md_dt <- data.table(
    dt1 = nrow(remove.from),
    dt2 = nrow(remove.with),
    merged_dt = nrow(merged_dt),
    duplicates = duplicate_count,
    nas = na_count,
    removed = nrow(remove.from) - nrow(merged_dt)
  )
  setnames(md_dt, c("dt1", "dt2"), c(dt1_name, dt2_name))
  
  md_dt <- kable(md_dt)
  
  vebprint(md_dt, text = "Anti-union summary:")
  
  mdwrite(
    config$files$post_seq_md,
    text = paste0("Removing data from **", dt1_name, "** using **", dt2_name, "**"),
    data = md_dt
  )
  
  return(merged_dt)
}

inner_union <- function(dt1, dt2, by, na.rm = TRUE, verbose = FALSE) {
  catn("Inner joining data tables")
  
  if (!is.data.table(dt1) || !is.data.table(dt2)) {
    stop("Both dt1 and dt2 must be data.tables")
  }
  
  dt1_name <- deparse(substitute(dt1))
  dt2_name <- deparse(substitute(dt2))
  
  vebcat("Merging", highcat(dt1_name), "and", highcat(dt2_name), veb = verbose)
  
  merged_dt <- dt1[dt2, nomatch = 0, on = by]
  
  ## Run NA and duplicate check
  na_count <- sum(is.na(merged_dt))
  
  if (any(is.na(merged_dt)) == TRUE) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    if (na.rm) merged_dt <- na.omit(merged_dt)
  }
  
  duplicate_count <- sum(duplicated(merged_dt))
  
  if (any(duplicated(merged_dt)) == TRUE) {
    merged_dt <- unique(merged_dt)
  }
  
  md_dt <- data.table(
    dt1 = nrow(dt1),
    dt2 = nrow(dt2),
    merged_dt = nrow(merged_dt),
    duplicates = na_count,
    nas = duplicate_count
  )
  setnames(md_dt, c("dt1", "dt2"), c(dt1_name, dt2_name))
  
  md_dt <- kable(md_dt)
  
  vebprint(md_dt, text = "Inner join summary:")
  
  mdwrite(
    config$files$post_seq_md,
    text = paste0("Inner joining **", dt1_name, "** and **", dt2_name, "**"),
    data = md_dt
  )
  
  return(merged_dt)
}

select_wfo_column <- function(dir.path, col.unique, col.select = NULL, col.combine = NULL, pattern = "*.csv", verbose = FALSE) {
  catn("Selecting WFO column")
  
  csv_files <- get_files(dir.path, include.files = pattern)
  
  dt_list <- list()
  
  for (file in csv_files) {
    dt <- fread(file)
    
    split_str <- str_split(file, "/")[[1]]
    parent_name <- split_str[length(split_str) - 2]
    child_name <- split_str[length(split_str) - 1]
    name <- paste0(parent_name, "_", child_name)
    
    vebcat("Selecting columns for", highcat(name), veb = verbose)
    
    # Remove only genus names
    dt <- dt[tolower(taxonRank) != "genus"]
    
    if (!is.null(col.select)) {
      vebcat("Selecting the", highcat(col.select), "column(s) and using", highcat(col.unique), "as the unique column.", veb = verbose)
      dt_sel <- dt[, c(col.select, col.unique), with = FALSE]
    } else {
      vebcat("Using the", col.unique, "as unique column.", veb = verbose)
      dt_sel <- dt[, ..col.unique, with = FALSE]
    }
    
    dt_uniq <- dt_sel[!duplicated(dt_sel[[col.unique]]), ]
    
    if (verbose) {
      n_orig <- nrow(dt_sel)
      n_uniq <- nrow(dt_uniq)
      catn("\nList:", highcat(name))
      cat(sprintf("%-10s | %s \n", "n_species", "unique(n_species)"))
      cat(highcat(sprintf("%-10d | %d \n", n_orig, n_uniq)))
      catn()
    }
    
    # Add this data.table to our list with the constructed name
    dt_list[[name]] <- dt_uniq
  }
  
  return(dt_list = dt_list)
}

fix_manual <- function(dts, manual.edited, column, verbose = FALSE) {
  catn("Fixing manual edits.")
  
  if (!file.exists(manual.edited)) {
    vebcat("Could not detect manually edited file", color = "fatalError")
    catn("File can be found at:", colcat("./outputs/setup/wrangle/manual-check-file.csv", color = "output"))
    catn("Help can be found in the GitHub README:", highcat("https://github.com/TorHenrikUlsted/arctic-horizon-scanning"))
    
    stop("Edit manual file and then rerun")
  }
  
  # Combine no-matches
  edited <- fread(manual.edited, select = c("acceptedName", "listOrigin"), encoding = "UTF-8")
  
  # Remove NAs and blank values in the newName df
  formatted <- edited[!is.na(acceptedName) & acceptedName != "" & acceptedName != "removed"]
  
  formatted[, acceptedName.ORIG := acceptedName]
  
  formatted[, acceptedName := { # Clean symbols & designations
    tmp <- clean_string(acceptedName, verbose)
    tmp <- clean_designations(tmp, config$species$standard_infraEpithets, verbose)
    tmp <- clean_symbols(tmp, config$species$standard_symbols, verbose)
    acceptedName = gsub("× ", "×", tmp)
  }]
  
  data.table::setnames(formatted, old = "acceptedName", new = "scientificName")
  
  vebprint(formatted, verbose, "Manually formatted data:")
  
  for (i in 1:nrow(formatted)) {
    # Get the scientificName and listOrigin from the current row
    scientificName <- formatted[i, scientificName]
    listOrigin <- as.character(formatted[i, listOrigin])
    
    vebcat("Appending", 
           highcat(paste(as.character(scientificName))), 
           "to", highcat(listOrigin), veb = verbose
          )
    
    # Check if scientificName or listOrigin is missing
    if (is.na(scientificName) | is.na(listOrigin)) {
      catn("Missing value in row ", i)
      next
    }
    
    # Check if listOrigin exists in dts
    if (!listOrigin %in% names(dts)) {
      vebcat("Invalid listOrigin in row ", i, ": ", listOrigin, color = "nonFatalError")
      next
    }
    
    # Create a new data table with the same columns as dts[[listOrigin]]
    new_row <- data.table(scientificName = scientificName)
    
    setnames(new_row, names(dts[[listOrigin]]))
    
    # Append the scientificName to the corresponding data table in dts
    dts[[listOrigin]] <- rbindlist(list(dts[[listOrigin]], new_row), fill = TRUE)
  }
  
  vebcat("the manually formatted synonym checks have been successfully added to correct data frames.", veb = verbose, color = "proSuccess")
  
  return(dts)
}

write_filter_fun <- function(file.out, spec.in, fun = NULL) {
  create_dir_if(dirname(file.out))
  
  if (!is.null(fun)) {
    result <- fun()
  } else {
    result <- spec.in
  }
  
  catn("Writing species to:", colcat(file.out, color = "output"))
  
  fwrite(result, file.out, bom = T)
  
  return(result)
}

# filter known and unknown using GBIF keys
filter_gbif_keys <- function(spec.dts, out.dirs, verbose = FALSE) {
  
  present_out <- paste0(out.dirs$unknown, "/present-usage-keys.csv")
  absent_out <- paste0(out.dirs$unknown, "/absent-usage-keys.csv")
  
  if (file.exists(present_out) & file.exists(absent_out)) return(fread(absent_out))
  
  mdwrite(
    config$files$post_seq_md,
    text = "2;Filtering using GBIF keys"
  )
  
  known_keys <- gbif_standardize(
    dt = spec.dts$known,
    out.file = paste0(out.dirs$known, "/standardized-sp-keys.csv"),
    verbose = verbose
  )
  
  data.table::setnames(known_keys, "usageKey.GBIF", "usageKey")
  known_keys <- known_keys[, .(usageKey)]
  
  unknown_keys <- gbif_standardize(
    dt = spec.dts$unknown,
    out.file = paste0(out.dirs$unknown, "/standardized-sp-keys.csv"),
    verbose = verbose
  )
  
  data.table::setnames(unknown_keys, "usageKey.GBIF", "usageKey")
  # Filter out usageKeys that will be included in speciesKey
  species_keys <- unknown_keys[tolower(rank.GBIF) == "species", unique(speciesKey.GBIF)]
  # Filter out lower-level taxa if their species is already present
  removed_keys <- unknown_keys[speciesKey.GBIF %in% species_keys & tolower(rank.GBIF) != "species"]
  
  if (nrow(removed_keys) > 0) {
    catn("Found", highcat(nrow(removed_keys)), "taxons below species where the species is already included")
    fwrite(removed_keys, paste0(out.dirs$unknown, "/dup_sp_keys.csv"))
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0("**", nrow(removed_keys), "** taxons below an already included species were removed")
    )
  }
  
  unknown_keys <- unknown_keys[!(speciesKey.GBIF %in% species_keys & tolower(rank.GBIF) != "species")]
  
  fwrite(unknown_keys, paste0(out.dirs$unknown, "/standardized-sp-keys.csv"))
  
  unknown_keys <- unknown_keys[, .(usageKey)]
  
  unknown_present <- write_filter_fun(
    file.out = present_out,
    spec.in = unknown_keys,
    fun = function() {
      # First merge to only get species from both dts
      unknown_present <- inner_union(unknown_keys, known_keys, by = "usageKey")
      
      return(unknown_present)
    }
  )
  
  unknown_absent <- write_filter_fun(
    file.out = absent_out,
    spec.in = unknown_keys,
    fun = function() {
      # Remove known present species from the unknown-absent species
      unknown_absent <- anti_union(unknown_keys, known_keys, "usageKey")
      
      return(unknown_absent)
    }
  )
  
  mdwrite(
    config$files$post_seq_md,
    text = paste0("Number of species keys for download: **", nrow(unknown_absent), "**")
  )
  
  return(unknown_absent)
}

get_prefix <- function(taxonRank) {
  switch(taxonRank,
         "SPECIES" = "",
         "SUBSPECIES" = "subsp.",
         "VARIETY" = "var.",
         "FORM" = "f.",
         ""
  )
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

# spec.occ can be either filepath or data frame
chunk_protocol <- function(
    spec.occ,
    spec.keys = NULL,
    chunk.name = "species",
    chunk.col,
    chunk.dir,
    chunk.size = 1e6,
    iterations = NULL,
    approach = FALSE,
    verbose = FALSE) {
  
  vebprint(head(spec.occ, 1), text = "spec.occ", veb = verbose)
  
 chunk_data(
   spec.occ = spec.occ, 
   chunk.name = chunk.name, 
   chunk.column = chunk.col, 
   chunk.dir = chunk.dir, 
   chunk.size = chunk.size, 
   iterations = iterations, 
   verbose = verbose
 )
 
  rename_chunks(
    spec.dir = file.path(chunk.dir, "species"),
    symbols = config$species$filename_symbols,
    designations = config$species$standard_infraEpithets,
    file.sep = config$species$file_separator,
    verbose = verbose
  )
  
  if (approach == "conservative" & !is.null(spec.keys)) { # this has to be fixed specifically for new conservative approach.
    clean_chunks(
      chunk.name = chunk.name,
      chunk.column = chunk.col,
      chunk.dir = chunk.dir,
      sp_w_keys = spec.keys,
      verbose = verbose
    )
  }
}
