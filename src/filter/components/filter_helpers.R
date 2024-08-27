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

  ## Run NA and distinct check
  if (any(is.na(merged_dt)) == T) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    if (na.rm) merged_dt <- na.omit(merged_dt)
  }

  if (any(duplicated(merged_dt)) == T) {
    vebcat("Some merged_dt species are duplicated.", color = "nonFatalError")
    merged_dt <- unique(merged_dt)
  }

  duplicate_count <- sum(duplicated(merged_dt))
  na_count <- sum(is.na(merged_dt))

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

  vebprint(md_dt, text = "Anti-join summary:")

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
  if (any(is.na(merged_dt)) == TRUE) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    if (na.rm) merged_dt <- na.omit(merged_dt)
  }

  if (any(duplicated(merged_dt)) == TRUE) {
    vebcat("Some merged_dt species are duplicated.", color = "nonFatalError")
    merged_dt <- unique(merged_dt)
  }

  duplicate_count <- sum(duplicated(merged_dt))
  na_count <- sum(is.na(merged_dt))

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

  vebprint(md_dt, text = "Anti-join summary:")

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
  if (any(is.na(merged_dt)) == TRUE) {
    vebcat("Some merged data table species are NA.", color = "nonFatalError")
    if (na.rm) merged_dt <- na.omit(merged_dt)
  }

  if (any(duplicated(merged_dt)) == TRUE) {
    vebcat("Some merged_dt species are duplicated.", color = "nonFatalError")
    merged_dt <- unique(merged_dt)
  }

  duplicate_count <- sum(duplicated(merged_dt))
  na_count <- sum(is.na(merged_dt))

  md_dt <- data.table(
    dt1 = nrow(dt1),
    dt2 = nrow(dt2),
    merged_dt = nrow(merged_dt),
    duplicates = na_count,
    nas = duplicate_count,
    removed = (nrow(dt1) + nrow(dt2)) - nrow(merged_dt)
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

  # make a list of data frames based on the different CSV files and also check for any "no matches" then add those to their own data frame.
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

fix_nomatches <- function(dts, nomatch.edited, column, verbose = FALSE) {
  catn("Fixing nomatches.")

  # Combine no-matches
  edited_nomatch <- fread(nomatch.edited, encoding = "UTF-8")

  # Remove NAs and blank values in the newName df
  formatted <- edited_nomatch[!is.na(acceptedName) & acceptedName != "", .(acceptedName, listOrigin)]

  vebprint(formatted, verbose, "Manually formatted data:")

  for (i in 1:nrow(formatted)) {
    # Get the scientificName and listOrigin from the current row
    scientificName <- formatted[i, acceptedName]
    listOrigin <- as.character(formatted[i, listOrigin])

    vebcat("Appending", highcat(as.character(scientificName)), "to", highcat(listOrigin), veb = verbose)

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

remove_authorship <- function(dt, verbose = FALSE) {
  vebcat("Removing Authorship from scientificName.", veb = verbose)

  signs <- c("×")

  dt[, cleanName := sapply(scientificName, function(x) {
    components <- strsplit(x, " ")[[1]]
    result <- c(components[1], components[2])

    if (components[2] %in% signs) {
      result <- c(result, components[3])
    }

    if (components[3] %in% config$species$infraEpithet_designations) {
      result <- c(result, components[3], components[4])
    }

    if (components[4] %in% config$species$infraEpithet_designations || components[4] %in% signs) {
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
    approach = FALSE,
    verbose = FALSE) {
  vebprint(head(spec.occ, 1), text = "spec.occ", veb = verbose)

  if (is.character(spec.occ)) {
    chunk_file(
      file.path = spec.occ,
      approach = approach,
      chunk.name = chunk.name,
      chunk.column = chunk.col,
      chunk.dir = chunk.dir,
      chunk.size = chunk.size,
      cores.max = cores.max,
      iterations = iterations,
      verbose = verbose
    )
  } else if ("data.frame" %in% class(spec.occ) || "data.table" %in% class(spec.occ)) {
    chunk_loaded_df(
      dt = spec.occ,
      approach = approach,
      chunk.name = chunk.name,
      chunk.column = chunk.col,
      chunk.dir = chunk.dir,
      verbose = verbose
    )
  } else {
    stop("Invalid input, either filepath or data frame/table")
  }

  if (approach == "conservative") { # this has to be fixed specifically for new conservative approach.
    clean_chunks(
      chunk.name = chunk.name,
      chunk.column = chunk.col,
      chunk.dir = chunk.dir,
      sp_w_keys = spec.keys,
      verbose = verbose
    )
  }
}
