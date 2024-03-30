##########################
#       Data.table       #
##########################

merge_and_sum <- function(dt1, dt2, sumCol, by, all = TRUE) {
  merged_dt <- merge(dt1, dt2, by = by, all = all)
  merged_dt[[sumCol]] <- rowSums(merged_dt[, c(paste0(sumCol, ".x"), paste0(sumCol, ".y"))], na.rm = TRUE)

  return(merged_dt)
}

find_min_data_col <- function(file_path, verbose = TRUE) {
  catn("Finding column with the least memory allocation needed.")
  # Read the first 10 rows of the file
  dt <- fread(file_path, nrows = 10)

  # Calculate the total character length of each column
  column_lengths <- sapply(dt, function(x) sum(nchar(as.character(x))))

  # Find the column name with the least character length
  least_data_column <- names(column_lengths)[which.min(column_lengths)]

  catn("Column with the least data:", highcat(least_data_column))

  return(least_data_column)
}

set_df_utf8 <- function(df) {
  for (name in names(df)[sapply(df, is.character)]) {
    df[[name]] <- enc2utf8(df[[name]])
  }

  return(df)
}

remove_infraEpithet <- function(spec) {
  for (d in infraEpithet_designations) {
    # Remove the designation and the name that follows it from the species name
    spec <- gsub(paste0("(\\s*\\(?\\s*(?i)", d, "\\.?\\s+)([^\\)]*)\\)?"), "", spec)
  }

  return(spec)
}

# Get the infraspecificEpithet and standardize it
extract_infraEpithet <- function(spec, verbose = FALSE) {
  res <- ""

  for (d in infraEpithet_designations) {
    pattern <- paste0("(?i)\\b", gsub("\\.", "", d), "\\.?\\s+(\\w+)\\b")

    if (grepl(pattern, spec)) {
      match <- regmatches(spec, regexpr(pattern, spec, perl = TRUE))

      if (length(match) > 0) {
        res <- match[[1]]
        break
      }
    }
  }

  vebprint(res, veb = verbose)

  # Standardize the infraspecific epithet
  if (res != "") {
    for (s in names(standard_infraEpithets)) {
      pattern <- paste0("(?i)\\b", gsub("\\.", "", s), "\\.?\\s+(\\w+)\\b")

      if (grepl(pattern, res)) {
        epithet <- tolower(gsub(pattern, "\\1", res, perl = TRUE))
        separator <- ifelse(grepl("\\.$", s), "", ". ")
        replacement <- paste0(standard_infraEpithets[[s]], separator, " ", epithet)
        res <- replacement

        break
      }
    }
  }

  vebprint(res, veb = verbose)

  return(res)
}

# Used to only standardize, can be used with sapply for certain columns in df
standardize_infraEpithet <- function(spec, verbose = FALSE) {
  res <- spec

  for (s in names(standard_infraEpithets)) {
    pattern <- paste0("(?i)\\b", gsub("\\.", "", s), "\\.?\\s+(\\w+)\\b")

    if (grepl(pattern, res)) {
      epithet <- tolower(gsub(pattern, "\\1", res, perl = TRUE))
      separator <- ifelse(grepl("\\.$", s), "", ". ")
      replacement <- paste0(standard_infraEpithets[[s]], separator, " ", epithet)
      res <- replacement

      break
    }
  }

  return(res)
}

select_wfo_column <- function(filepath, col.unique, col.select = NULL, col.combine = NULL, pattern = "*.csv", verbose = F) {
  
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
    catn("\nList:", highcat(sub("-wfo-one.csv$", "", basename(file))))
    cat(sprintf("%-10s | %s \n", "n_species", "unique(n_species)"))
    cat(highcat(sprintf("%-10d | %d \n", n_orig, n_uniq)))
    catn()
    
    return(df_uniq)
  })
  
  names(df_list) <- sub("-wfo-one.csv$", "", basename(csv_files))
  names(df_list) <- gsub("-", "_", names(df_list))
  
  
  return(df_list = df_list)
}

##########################
#         Terra          #
##########################

extract_ext_to_dt <- function(raster, value = "value", cells = TRUE) {
  rast_extr <- terra::extract(raster, ext(raster), cells = cells)
  rast_dt <- as.data.table(rast_extr)
  names(rast_dt) <- c("cell", value)

  return(rast_dt)
}

convert_spatial_dt <- function(spatial, verbose = FALSE) {
  dt <- as.data.frame(spatial)
  dt <- as.data.table(dt)

  return(dt)
}

reproject_region <- function(region, projection, line_issue = F, show_plot = F, verbose = T) {
  region_name <- strsplit(deparse(substitute(region)), "\\$")[[1]][[2]]
  catn("Reprojecting", region_name)

  if (line_issue == T) {
    catn("Attempting to fix line issues.")

    catn("Getting extents.")
    ext_east <- terra::ext(ext(region)$xmin, 0, ext(region)$ymin, ext(region)$ymax)
    ext_west <- terra::ext(0.00001, ext(region)$xmax, ext(region)$ymin, ext(region)$ymax)

    catn("Cropping in half.")
    vect_east <- terra::crop(region, ext_east)
    vect_west <- terra::crop(region, ext_west)

    catn("Reprojecting to longlat.")
    proj_east <- terra::project(vect_east, longlat_crs)
    catn("Plotting left side.")
    plot(proj_east)
    proj_west <- terra::project(vect_west, longlat_crs)
    catn("plotting right side.")
    plot(proj_west)

    region_longlat <- rbind(proj_west, proj_east)

    if (show_plot == T) plot(region)

    return(region_longlat)
  }


  if (projection == "longlat") {
    catn("Choosing", highcat("longlat"), "coordinate system.")

    if (grepl("+proj=longlat", crs(region, proj = T), fixed = TRUE) == FALSE) {
      catn("Is not longlat, will be reprojected.")
    }

    print(crs(longlat_crs, proj = T))
    prj <- longlat_crs
  } else if (projection == "laea") {
    catn("Choosing", highcat("polar"), "coordinate system. \n")
    if (grepl("+proj=laea", crs(region, proj = T), fixed = TRUE) == FALSE) {
      vebcat("Is not polar.", color = "nonFatalError")
    }
    prj <- laea_crs
  } else {
    stop("You can only choose projection 'longlat' or 'laea'.")
  }


  if (!isTRUE(identical(crs(region), prj))) {
    vebcat("Original CRS not identical to current CSR.", veb = verbose)
    vebcat("Reprojecting", highcat(region_name), "to: ", crs(prj, proj = T), veb = verbose)

    reproj_region <- terra::project(region, prj)

    if (!isTRUE(identical(crs(reproj_region), prj))) {
      vebcat("Reprojection completed successfully", color = "funSuccess")
    } else {
      vebcat("Reprojection failed.", color = "nonFatalError")
    }
  } else {
    vebcat("Original CRS identical to current CSR.", veb = verbose)
    vebcat(region_name, "reprojected successfully.", color = "funSuccess")
  }

  return(reproj_region)
}

fix_shape <- function(shape, verbose = FALSE) {
  if (any(!is.valid(shape))) {
    catn("Some geoms of", substitute(deparse(shape)), "are invalid.")
    catn("Attempting to fix.")
    valid_shape <- makeValid(shape)
    if (any(!is.valid(valid_shape))) {
      stop("Failed to fix invalid geoms.")
    } else {
      vebcat("Successfully made all geoms valid.", color = "proSuccess")
    }
    return(valid_shape)
  } else {
    catn("All", substitute(deparse(shape)), "geoms are valid.")

    return(shape)
  }
}

order_by_apg <- function(input, by, verbose = FALSE) {
  # Load apg4
  apg <- fread("./resources/taxon/apg4/apg4.txt")

  apg_rank <- apg[taxonRank == by]

  apg_vect <- apg_rank$scientificName

  non_apg <- setdiff(input, apg_vect)

  vebcat("non_apg:", veb = verbose)
  vebprint(non_apg, veb = verbose)

  input_apg <- input[!input %in% non_apg]

  vebcat("input without non_apgs:", veb = verbose)
  vebprint(input_apg, veb = verbose)

  vebcat("Ordering apgs", veb = verbose)
  if (is.vector(input)) {
    ordered_apg <- apg_vect[apg_vect %in% input_apg]

    # Add the others
    non_apg_order <- c("Lycopodiales", "Selaginellales", "Equisetales", "Osmundales", "Salviniales", "Cyatheales", "Polypodiales", "Ephedrales", "Pinales")

    ordered <- c(non_apg_order, ordered_apg)
    print(ordered)
  } else {
    setorderv(input, cols = apg_rank[[by]])
  }

  return(ordered)
}

find_peaks <- function(data, prominence = 0.1) {
  # Identify local maxima using diff
  peaks <- which(diff(data) > 0 & diff(c(data, 0)) < 0)

  # Filter based on prominence (difference from neighbors)
  filtered_peaks <- peaks[data[peaks] - data[peaks - 1] >= prominence & data[peaks] - data[peaks + 1] >= prominence]

  return(filtered_peaks)
}

##########################
#       Parallel         #
##########################
load_sp_rast <- function(spec.filename) {
  sp_name <- basename(dirname(spec.filename))

  sp_rast <- terra::rast(spec.filename)
  names(sp_rast) <- sp_name

  return(sp_rast)
}

##########################
#      File system       #
##########################

source_all <- function(dir) {
  # Get a list of all .R files in the directory and its subdirectories
  r_files <- list.files(path = dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)

  # Source each file
  lapply(r_files, source)

  cat(length(r_files), "scripts sourced from", dir, "and its subdirectories.\n")
}
