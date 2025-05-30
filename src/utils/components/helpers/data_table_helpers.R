#------------------------#
####    Data.table    ####
#------------------------#

merge_and_sum <- function(dt1, dt2, sumCol, by, all = TRUE) {
  merged_dt <- merge(dt1, dt2, by = by, all = all)
  merged_dt[[sumCol]] <- rowSums(merged_dt[, c(paste0(sumCol, ".x"), paste0(sumCol, ".y"))], na.rm = TRUE)
  
  return(merged_dt)
}

find_min_data_col <- function(dt, count.rows = 100, verbose = TRUE) { # count.rows should be sample
  catn("Finding column with the least memory allocation needed.")
  
  byte_size <- list(
    logical = 1,
    integer = 4,
    integer64 = 8,
    numeric = 8,
    double = 8,
    POSIXct = 8,
    POSIXt = 8,
    IDate = 4,
    Date = 4,
    character = function(string) {
      1 + nchar(string)
    }
  )
  
  # Read the first count.rows number of rows of the file
  if (!(is.data.table(dt) || is.data.frame(dt)) && is.character(dt)) {
    dt <- fread(dt, nrows = count.rows)
  } else {
    return(vebcat("The input is not  a filepath nor a data frame or data table.", color = "fatalError"))
  }
  
  # Calculate the total byte size of each column
  column_sizes <- sapply(dt, function(x) {
    for (class_name in names(byte_size)) {
      if (inherits(x, class_name)) {
        if (class_name == "character") {
          return(sum(byte_size[[class_name]](x)))
        } else {
          return(length(x) * byte_size[[class_name]])
        }
      }
    }
    warning("Warning: Class '", paste(class(x), collapse = ", "), "' is not in the byte_size list. Skipping this column.")
    return(NA) # Return NA for classes not in the list
  })
  
  column_sizes <- column_sizes[!is.na(column_sizes)]
  
  # Find the column name with the least data
  least_data_column <- names(column_sizes)[which.min(column_sizes)]
  
  catn("Column with the least data:", highcat(least_data_column))
  
  rm(dt, column_sizes)
  invisible(gc())
  
  return(least_data_column)
}

set_df_utf8 <- function(df) {
  for (name in names(df)[sapply(dt, is.character)]) {
    df[[name]] <- enc2utf8(df[[name]])
  }
  
  return(df)
}

# Used to only standardize, can be used with sapply for certain columns in df
standardize_infraEpithet <- function(spec, verbose = FALSE) {
  res <- spec
  
  for (s in names(config$species$standard_infraEpithets)) {
    pattern <- paste0("(?i)\\b", gsub("\\.", "", s), "\\.?\\s+(\\w+)\\b")
    
    if (grepl(pattern, res)) {
      epithet <- tolower(gsub(pattern, "\\1", res, perl = TRUE))
      separator <- ifelse(grepl("\\.$", s), "", ". ")
      replacement <- paste0(config$species$standard_infraEpithets[[s]], separator, " ", epithet)
      res <- replacement
      
      break
    }
  }
  
  return(res)
}

remove_designations <- function(spec, verbose = FALSE) {
  # Remove ignored designations
  for (d in config$species$ignored_designations) {
    pattern <- paste0("\\b", d, "\\b(?:\\s*\\([^)]+\\))?")
    
    spec <- gsub(pattern, "", spec)
  }
  
  return(trimws(spec))
}

remove_infraEpithet <- function(spec, verbose = FALSE) {
  # Remove ignored designations
  spec <- remove_designations(spec = spec, verbose = verbose)
  
  for (d in config$species$infraEpithet_designations) {
    # Remove the designation and the name that follows it from the species name
    spec <- gsub(paste0("(\\s*\\(?\\s*(?i)", d, "\\.?\\s+)([^\\)]*)\\)?"), "", spec)
  }
  
  return(spec)
}

# Get the infraspecificEpithet and standardize it
extract_infraEpithet <- function(spec, verbose = FALSE) {
  res <- ""
  
  for (d in config$species$infraEpithet_designations) {
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
    res <- standardize_infraEpithet(res, verbose = verbose)
  }
  
  vebprint(res, veb = verbose)
  
  return(res)
}

reverse_list <- function(list) {
  reversed <- setNames(as.list(names(list)), unlist(list))
  return(reversed)
}

clean_string <- function(x, verbose = FALSE) {
  x <- gsub("–", "-", x) # dashes to hyphens
  x <- gsub("\u00A0", " ", x) # Replace non-breaking spaces with regular spaces
  x <- gsub("([A-Z]),", "\\1.", x) # If one capital letter then ",", change to ".", then remove double "."
  x <- gsub('"', "", x) # Remove " in a string
  x <- gsub("\\.\\.", ".", x) # remove double "." with only one
  # Remove standalone "x" where the next word is longer than 1 character and ends with a period
  x <- gsub("\\s+x\\s+(?=\\S{2,}\\.)\\s*", " ", x, perl = TRUE)
  
  x <- gsub("\\s+", " ", x) # Remove double spaces
  x <- trimws(x) # trim trailing spaces
  
  # vebcat("Cleaned string:", x, veb = verbose)
  # vebcat("ASCII values:", paste(sapply(strsplit(x, "")[[1]], function(ch) as.integer(charToRaw(ch))), collapse = " "), veb = verbose)
  
  return(x)
}

create_symbol_pattern <- function(symbol, verbose = FALSE) {
  escaped_symbol <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", symbol)
  
  patterns <- if (tolower(symbol) == "x") {
    c(
      paste0("(?<=\\s|^)[", escaped_symbol, toupper(escaped_symbol), "](?=\\s|$)"), # standalone
      paste0("(?<=\\s)[", escaped_symbol, toupper(escaped_symbol), "](?=\\s)") # surrounded by spaces
    )
  } else if (symbol == "×") {
    escaped_symbol # anywhere in the string
  } else {
    c(
      paste0("(?<=\\s|^)", escaped_symbol, "(?=\\s|$)"), # standalone
      paste0("(?<=\\s)", escaped_symbol, "(?=\\s)") # surrounded by spaces
    )
  }
  
  if (is.character(patterns) && length(patterns) == 1) {
    patterns <- list(patterns) # Ensure patterns is always a list
  }
  
  return(patterns)
}

clean_symbols <- function(x, symbols, verbose = FALSE) {
  for (symbol in names(symbols)) {
    replacement <- symbols[[symbol]]
    
    patterns <- create_symbol_pattern(symbol, verbose)
    
    for (pattern in patterns) {
      x <- gsub(pattern, paste0(" ", replacement, " "), x, perl = TRUE)
      # vebcat(paste("Replaced", pattern, "with", replacement), veb = verbose)
    }
  }
  
  x <- trimws(gsub("\\s+", " ", x)) # Clean up any extra spaces
  return(x)
}

create_designation_pattern <- function(designation, verbose = FALSE) {
  # Remove trailing period if present
  base_designation <- sub("\\.$", "", designation)
  
  # Escape special characters in the designation
  escaped_designation <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", base_designation)
  
  patterns <- c(
    paste0("(?<=\\s|^)", escaped_designation, "(?=\\s|$)"), # without period
    paste0("(?<=\\s|^)", escaped_designation, "\\.?(?=\\s|$)"), # with optional period
    paste0("(?<=\\s|^)", tolower(escaped_designation), "(?=\\s|$)"), # lowercase without period
    paste0("(?<=\\s|^)", tolower(escaped_designation), "\\.?(?=\\s|$)") # lowercase with optional period
  )
  
  return(patterns)
}

clean_designations <- function(x, designations, verbose = FALSE) {
  x <- gsub("\\.", ". ", x) # Add space after every "."
  
  for (designation in names(designations)) {
    # Create patterns for different cases
    patterns <- create_designation_pattern(designation, verbose)
    
    for (pattern in patterns) {
      # Apply the pattern and replacement
      x <- gsub(pattern, designations[[designation]], x, perl = TRUE)
    }
  }
  x <- gsub("\\. ", ".", x) # Remove space after every "."
  
  # Remove any double periods that might have been introduced
  x <- gsub("\\.\\.", ".", x)
  
  return(x)
}

uniq_list_pattern <- function(list, verbose = FALSE) {
  # Get unique standardized designations and symbols
  unique_items <- unique(unlist(list))
  # Escape designations and make pattern
  escaped_items <- gsub("([.|()\\^{}+$*?])", "\\\\\\1", unique_items)
  pattern <- paste(escaped_items, collapse = "|")
  
  return(pattern)
}

identify_structure <- function(x, symbol.pattern, designation.pattern, verbose = FALSE) {
  # Split the name, keeping parentheses together
  x <- gsub("\\(\\s", "(", x) # Remove space immediately after opening parenthesis
  x <- gsub("\\s\\)", ")", x) # Remove space immediately before closing parenthesis
  
  # Split the name, keeping parentheses together
  parts <- str_match_all(x, "<<[^>]+>>|\\([^()]*\\)|\\S+")[[1]]
  parts <- parts[parts != ""] # Remove any empty parts
  
  # Find the index of the designation, if it exists
  designation_indices <- which(sapply(parts, function(part) grepl(paste0("^(", designation.pattern, ")$"), part, perl = TRUE)))
  # Find the indices of the symbols, if they exist
  symbol_indices <- which(sapply(parts, function(part) grepl(paste0("^(", symbol.pattern, ")$"), part, perl = TRUE)))
  
  # vebprint(parts, verbose, text = "Parts:")
  # vebcat("symbol_indices:", symbol_indices, veb = verbose)
  # vebcat("symbol count:", length(symbol_indices), veb = verbose)
  # vebcat("designation_indices:", designation_indices, veb = verbose)
  # vebcat("designation_count:", length(designation_indices), veb = verbose)
  
  # Check if the designation is followed by "NA" at the end
  if (length(designation_indices) > 0) {
    to_remove <- c()
    for (i in 1:length(designation_indices)) {
      d <- designation_indices[i]
      
      # Check if this is the last part or if the next part is NA
      if (length(parts) == d || is.na(parts[d + 1])) {
        to_remove <- c(to_remove, i)
      }
      
      # Check if the word before the designation ends on a "." and is "f."
      if ((grepl("\\.$", parts[(d - 1)]) | grepl("Baker", parts[(d - 1)])) &
          grepl("^f\\.$", parts[d])) {
        to_remove <- c(to_remove, i)
      }
      
      # Check if the next index is also a designation
      if (i < length(designation_indices) && designation_indices[(i + 1)] == d + 1) {
        to_remove <- c(to_remove, i)
      }
      
      # Check if the word after the "f." is a single letter and not identical to genus
      if (grepl("^f\\.$", parts[d]) && grepl("^[A-Za-z]$", parts[(d + 1)]) && parts[(d + 1)] != substr(parts[1], 1, 1)) {
        to_remove <- c(to_remove, i)
      }
    }
    
    # Remove the flagged indices
    if (length(to_remove) > 0) {
      if (length(designation_indices) > 1) {
        designation_indices <- designation_indices[-to_remove]
      } else {
        # If there's only one designation and it's flagged for removal
        designation_indices <- NULL
      }
    }
    
    # If all designations were removed, set to NULL
    if (length(designation_indices) == 0) {
      designation_indices <- NULL
    }
  }
  
  # Remove symbols right before or after a designation
  if ((length(designation_indices) > 0 & length(symbol_indices) > 0) &
      (any(symbol_indices + 1 == designation_indices) || any(symbol_indices - 1 == designation_indices))
  ) {
    symbols_to_remove <- c()
    for (i in seq_along(designation_indices)) {
      for (j in seq_along(symbol_indices)) {
        if (symbol_indices[j] + 1 == designation_indices[i] || symbol_indices[j] - 1 == designation_indices[i]) {
          symbols_to_remove <- c(symbols_to_remove, symbol_indices[j])
        }
      }
    }
    
    if (length(symbols_to_remove) > 0) {
      parts <- parts[-symbols_to_remove]
      symbol_indices <- symbol_indices[!symbol_indices %in% symbols_to_remove]
    }
  }
  
  # CASES
  if (any(grepl("'([A-Z][^']+)'", parts)) & length(designation_indices) == 0) {
    structure <- fcase(
      length(symbol_indices) == 0 & length(designation_indices) == 0, "cultivar",
      length(symbol_indices) > 0 & length(designation_indices) == 0, "hybridCultivar",
      default = "unknownCultivar"
    )
  } else if (length(designation_indices) > 0 & length(symbol_indices) == 0) {
    structure <- fcase(
      length(designation_indices) == 1, "infraspecificTaxon",
      length(designation_indices) > 1, "MultiInfraspecificTaxon",
      default = "unknownInfraspecific"
    )
  } else if (length(designation_indices) == 0 & length(symbol_indices) > 0) {
    structure <- fcase(
      symbol_indices[1] == 1 || symbol_indices[1] > 1 &&
        symbol_indices[1] < length(parts) &&
        substr(parts[1], 1, 1) == toupper(substr(parts[1], 1, 1)) &&
        substr(parts[symbol_indices[1] + 1], 1, 1) == toupper(substr(parts[symbol_indices[1] + 1], 1, 1)) &&
        !grepl("^[A-Z]\\.$", parts[symbol_indices[1] + 1]) && # Not an abbreviated genus
        parts[1] != parts[symbol_indices[1] + 1], "interGenericHybrid",
      symbol_indices[1] == 2, "intraGenericHybrid",
      length(parts) >= 4 &&
        length(symbol_indices) == 1 &&
        symbol_indices[1] > 2 &&
        (parts[symbol_indices[1] + 1] == parts[1] ||
           grepl("^[A-Z]\\.$", parts[symbol_indices[1] + 1]) ||
           substr(parts[symbol_indices[1] + 1], 1, 1) == tolower(substr(parts[symbol_indices[1] + 1], 1, 1)) ||
           parts[symbol_indices[1] + 1] == tolower(parts[symbol_indices[1] + 1])), "interSpecificHybrid",
      length(parts) >= 4 && length(symbol_indices) > 1 && (symbol_indices[1] == 3 || symbol_indices[1] > 3), "MultiIntraSpecificHybrid",
      default = "unknownHybrid"
    )
  } else if (length(designation_indices) == 0 & length(symbol_indices) == 0) {
    structure <- "species"
  } else if (length(designation_indices) > 0 & length(symbol_indices) > 0) {
    structure <- fcase(
      any(grepl(symbol.pattern, parts[1])), "interGenericHybridTaxon",
      any(grepl(symbol.pattern, parts[2])) | any(grepl(symbol.pattern, parts[3])), "interSpecificHybridTaxon",
      default = "unknownInfraHybrid"
    )
  } else {
    structure <- "unknown"
  }
  
  # Initialize all fields with NA
  part_id <- list(
    genus = NA_character_,
    specificEpithet = NA_character_,
    infraspecificRank = NA_character_,
    hybrid = NA_character_,
    infraspecificEpithet = NA_character_,
    cleanName = NA_character_,
    other = NA_character_,
    fullName = NA_character_,
    extra = NA_character_,
    structure = structure
  )
  
  # Assign names based on the structure
  used_parts <- 0
  switch(structure,
         "species" = {
           part_id$genus <- parts[1]
           if (length(parts) >= 2 && grepl("\\.$", parts[2])) {
             part_id$other <- parts[2]
             parts <- parts[-2]
           }
           part_id$specificEpithet <- parts[2]
           part_id$cleanName <- paste(parts[1:2], collapse = " ")
           used_parts <- 2
         },
         "infraspecificTaxon" = {
           if (grepl("^[A-Za-z]\\.$", parts[(designation_indices + 1)])) {
             parts <- parts[-(designation_indices + 1)] # Remove the part after designation
           }
           part_id$genus <- parts[1]
           part_id$specificEpithet <- parts[2]
           part_id$infraspecificRank <- parts[designation_indices]
           infra_index <- 1
           if (designation_indices < (length(parts) - 1) &&
               grepl("^cv\\.|cvgr\\.$", parts[designation_indices]) &&
               grepl("^'", parts[designation_indices + 1]) &&
               !grepl("'$", parts[designation_indices + 1])) {
             infra_index <- 2
           }
           part_id$infraspecificEpithet <- parts[designation_indices + infra_index]
           part_id$cleanName <- paste(parts[c(1:2, designation_indices:(designation_indices + infra_index))], collapse = " ")
           used_parts <- designation_indices + infra_index
         },
         "MultiInfraspecificTaxon" = {
           part_id$extra <- paste(parts[c(1:2, designation_indices[2]:length(parts))], collapse = " ")
           parts <- parts[-(designation_indices[2]:length(parts))]
           
           part_id$genus <- parts[1]
           part_id$specificEpithet <- parts[2]
           part_id$infraspecificRank <- parts[designation_indices[1]]
           part_id$infraspecificEpithet <- parts[designation_indices[1] + 1]
           part_id$cleanName <- paste(parts[c(1:2, designation_indices[1]:(designation_indices[1] + 1))], collapse = " ")
           
           used_parts <- designation_indices[1] + 1
         },
         "interGenericHybrid" = {
           if (grepl("^[A-Z]", parts[1]) & substr(parts[1], 1, 1) != substr(parts[symbol_indices[1] + 1], 1, 1)) {
             if (symbol_indices != 3) {
               parts <- parts[-(3:(symbol_indices - 1))] # remove up to symbol
               symbol_indices <- symbol_indices - length(3:(symbol_indices - 1)) # set symbol to new position
             }
             part_id$genus <- parts[1]
             part_id$specificEpithet <- paste(parts[2:5], collapse = " ")
             part_id$hybrid <- parts[3]
             part_id$cleanName <- paste(parts[1:5], collapse = " ")
             used_parts <- 5
           } else {
             part_id$hybrid <- parts[1]
             part_id$genus <- paste(parts[1:2], collapse = " ")
             part_id$specificEpithet <- parts[3]
             part_id$cleanName <- paste(parts[1:3], collapse = " ")
             used_parts <- 3
           }
         },
         "intraGenericHybrid" = {
           part_id$genus <- parts[1]
           part_id$hybrid <- parts[2]
           part_id$specificEpithet <- paste(parts[2:3], collapse = " ")
           part_id$cleanName <- paste(parts[1:3], collapse = " ")
           used_parts <- 3
         },
         "interSpecificHybrid" = {
           if (symbol_indices != 3) {
             parts <- parts[-(3:(symbol_indices - 1))] # remove up to symbol
             symbol_indices <- symbol_indices - length(3:(symbol_indices - 1)) # set symbol to new position
           }
           if ((grepl("^[A-Z]\\.$", parts[symbol_indices[1] + 1]) || parts[symbol_indices[1] + 1] == parts[1]) && substr(parts[symbol_indices[1] + 1], 1, 1) == substr(parts[1], 1, 1)) {
             parts <- parts[-(symbol_indices[1] + 1)]
           } # Remove uppercase letter with "." after if identical to first uppercase letter of genus
           part_id$genus <- parts[1]
           part_id$specificEpithet <- paste(parts[2:4], collapse = " ")
           part_id$hybrid <- parts[3]
           part_id$cleanName <- paste(parts[1:4], collapse = " ")
           used_parts <- 4
         },
         "MultiIntraSpecificHybrid" = {
           part_id$genus <- parts[1]
           part_id$specificEpithet <- paste(parts[2:6], collapse = " ")
           part_id$hybrid <- parts[3]
           part_id$cleanName <- paste(parts[1:6], collapse = " ")
           used_parts <- 6
         },
         "interGenericHybridTaxon" = {
           part_id$hybrid <- parts[1]
           part_id$genus <- paste(parts[1:2], collapse = " ")
           part_id$specificEpithet <- parts[3]
           part_id$infraspecificRank <- parts[4]
           part_id$infraspecificEpithet <- parts[5]
           part_id$cleanName <- paste(parts[1:5], collapse = " ")
           used_parts <- 5
         },
         "interSpecificHybridTaxon" = {
           if (symbol_indices == 2) {
             part_id$genus <- parts[1]
             part_id$hybrid <- parts[2]
             part_id$specificEpithet <- paste(parts[2:3], collapse = " ")
             if (designation_indices != 4) {
               parts <- parts[-(4:(designation_indices - 1))]
               designation_indices <- designation_indices - length(4:(designation_indices - 1))
             }
             part_id$infraspecificRank <- parts[4]
             part_id$infraspecificEpithet <- parts[5]
             part_id$cleanName <- paste(parts[1:5], collapse = " ")
             used_parts <- 5
           } else {
             part_id$genus <- parts[1]
             part_id$specificEpithet <- paste(parts[2:4], collapse = " ")
             part_id$hybrid <- parts[3]
             if (designation_indices != 5) {
               parts <- parts[-(5:(designation_indices - 1))]
               designation_indices <- designation_indices - length(5:(designation_indices - 1))
             }
             part_id$infraspecificRank <- parts[5]
             part_id$infraspecificEpithet <- parts[6]
             part_id$cleanName <- paste(parts[1:6], collapse = " ")
             used_parts <- 6
           }
         },
         "cultivar" = {
           if (grepl("\\(.*\\)", parts[2])) parts <- parts[-2]
           cultivar_index <- which(grepl("'", parts[1:length(parts)]))
           if (cultivar_index > 3) {
             parts <- parts[-(3:(cultivar_index - 1))]
             cultivar_index <- cultivar_index - length(3:(cultivar_index - 1))
           }
           part_id$genus <- parts[1]
           if (cultivar_index == 2) {
             part_id$specificEpithet <- parts[2]
             part_id$cleanName <- paste(parts[1:2], collapse = " ")
             used_parts <- 2
           } else {
             part_id$specificEpithet <- paste(parts[2:3], collapse = " ")
             part_id$cleanName <- paste(parts[1:3], collapse = " ")
             used_parts <- 3
           }
         },
         "hybridCultivar" = {
           part_id$genus <- parts[1]
           part_id$hybrid <- parts[2]
           part_id$specificEpithet <- parts[3]
           part_id$cleanName <- paste(parts[1:3], collapse = " ")
           used_parts <- 3
         },
         {
           # For complex cases, assign parts to named fields as much as possible
           part_id$genus <- parts[1]
           if (length(parts) > 1) part_id$specificEpithet <- parts[2]
           if (length(parts) > 2) part_id$hybrid <- parts[3]
           if (length(parts) > 3) part_id$infraspecificEpithet <- parts[4]
           part_id$cleanName <- paste(parts[1:min(4, length(parts))], collapse = " ")
           used_parts <- min(4, length(parts))
         }
  )
  
  # Add other parts
  if (used_parts < length(parts)) {
    other_parts <- parts[(used_parts + 1):length(parts)]
    part_id$other <- paste(other_parts, collapse = " ")
  } else if (!is.na(part_id$other) && part_id$other == "") {
    part_id$other <- NA_character_
  }
  
  # Clean parentheses with "=" inside them and remove "hybrid complex"
  part_id$other <- gsub("\\s*\\([^()]*=.*?\\)|\\s*hybrid complex", "", part_id$other)
  
  # Convert empty strings to NA and ensure character type
  part_id$other <- as.character(part_id$other) # Ensure character type
  part_id$other[trimws(part_id$other) == ""] <- NA_character_
  
  # Remove any NA values from cleanName
  part_id$cleanName <- gsub("\\s+", " ", gsub("NA", "", part_id$cleanName))
  part_id$cleanName <- trimws(part_id$cleanName)
  
  if (is.null(part_id$infraspecificRank) || length(part_id$infraspecificRank) == 0) {
    part_id$infraspecificRank <- NA_character_
  }
  
  part_id$fullName <- ifelse(!is.na(part_id$other),
                             paste(part_id$cleanName, part_id$other),
                             part_id$cleanName
  )
  
  return(part_id)
}

clean_spec_name <- function(x, symbols, designations, verbose = FALSE) {
  symbol_pattern <- uniq_list_pattern(symbols, verbose)
  designation_pattern <- uniq_list_pattern(designations, verbose)
  
  
  # Split designations into separate parts
  for (designation in designation_pattern) {
    x <- gsub(paste0("(\\s|^)(", designation, ")(\\s|$)"), "\\1\\2 ", x, perl = TRUE)
  }
  
  parts_id <- identify_structure(
    x = x,
    symbol_pattern,
    designation_pattern,
    verbose = verbose
  )
  
  return(parts_id)
}

# Combine the a row for all specified columns as strings. If custom.col and custom.list are used, the custom.col will be swapped out with the data from the custom.list.
combine_columns_dt <- function(..., dt, column.name = "combined", custom.col = NULL, custom.list = NULL, sep = " ", verbose = FALSE) {
  catn("Combining columns.")
  
  columns <- lapply(substitute(list(...))[-1], as.character)
  
  vebprint(columns, verbose, "Columns:")
  
  combined_dt <- copy(dt)
  
  vebprint(combined_dt, verbose, "Input data table:")
  
  # this will replace taxonRanks with the config taxonRanks
  if (!is.null(custom.col) && !is.null(custom.list)) {
    vebcat("custom.col and custom.list are not null.", veb = verbose, color = "proSuccess")
    
    if (custom.col %in% names(combined_dt)) {
      vebcat("custom.col is in the config.", veb = verbose, color = "proSuccess")
      combined_dt[, customList := custom.list[get(custom.col)]]
      
      pos <- match(custom.col, columns)
      
      vebprint(pos, verbose, "custom.col position:")
      
      columns[[pos]] <- "customList"
      
      vebprint(unique(combined_dt$customList), verbose, "Unique Custom list items:")
    }
  } else {
    vebcat("To use custom.column you need to specify the column name as a string.", color = "fatalError")
    vebcat("To use custom.list you need to specify the list as a list object.", color = "fatalError")
    stop("Both custom.col and custom.list have to be not null if they are to be used.")
  }
  
  combined_dt[, (column.name) := apply(
    .SD,
    1,
    function(x) paste(na.omit(x), collapse = sep)
  ),
  .SDcols = unlist(columns)
  ]
  
  combined_dt[, (column.name) := trimws(combined_dt[[column.name]])]
  
  vebprint(unique(combined_dt[[column.name]]), verbose, "Unique combinations:")
  
  return(combined_dt)
}

clean_spec_filename <- function(x) {
  
  base <- basename(x)
  
  clean_name <- tools::file_path_sans_ext(base)
  
  clean_name <- gsub(config$species$file_separator, " ", clean_name)
  
  return(clean_name)
}

get_spec_taxons <- function(spec) { 
  if (length(spec) == 1) {
    vebcat("Checking for one species")
    fun <- "name_backbone"
  } else if (length(spec) > 1) {
    vebcat("Using checklist")
    fun <- "name_backbone_checklist"
  } else if (nrow(spec) > 1) {
    vebcat("Using checklist")
    fun <- "name_backbone_checklist"
  } else {
    vebcat("Checking for one species")
    fun <- "name_backbone"
  }
  
  tryCatch({
    result <- do.call(fun, list(name = spec))
    if (any(is.na(result$order))) {
      vebcat("Some species returned NA:", color = "nonFatalError")
      index <- which(is.na(result$order))
      vebprint(spec[index])
    }
    return(result)
  }, error = function(e) {
    message("First attempt failed. Retrying...")
    Sys.sleep(0.5)  
    tryCatch({
      result <- do.call(fun, list(name = spec))
      if (any(is.na(result$order))) {
        vebcat("Some species returned NA:", color = "nonFatalError")
        index <- which(is.na(result$order))
        vebprint(spec[index])
      }
      return(result)
    }, error = function(e) {
      message("Second attempt also failed. Error: ", e$message)
      return(NULL)
    })
  })
}

get_spec_group <- function(spec) {
  res <- get_spec_taxons(spec)
  
  if (res$order %in% config$species$angiosperms) {
    result <- "angiosperms"
  } else if (res$order %in% config$species$gymnosperms) {
    result <- "gymnosperms"
  } else if (res$order %in% config$species$pteridophytes) {
    result <- "pteridophytes"
  } else {
    result <- "unknowns"
  }
  
  return(result)
}

get_spec_group_dt <- function(spec, spec.col = NULL, verbose = FALSE) {
  dt <- copy(spec)
  
  if (is.null(spec.col)) {
    spec.col <- names(dt)[1]
    vebcat("Assuming species names is in the first column in the data table, found:", highcat(spec.col))
  }
  
  
  dt[, group := fcase(
    is.na(order), "unknowns",
    order %in% config$species$angiosperms, "angiosperms",
    order %in% config$species$gymnosperms, "gymnosperms",
    order %in% config$species$pteridophytes, "pteridophytes",
    default = "others"
  )]
  
  vebprint(dt, verbose, "After adding groups:")
  
  unknown_species <- dt[group %in% c("unknowns", "others"), get(spec.col)]
  if (length(unknown_species) > 0) {
    vebcat("Order not found for species:", color = "nonFatalError")
    vebprint(unknown_species)
  }
  
  return(dt)
}

get_order_group <- function(dt, spec.col = NULL, verbose = FALSE) {
  dt_res <- get_spec_group_dt(dt, spec.col)
  
  # Concatenate the vectors in the desired order
  all_orders <- c(config$species$pteridophytes, config$species$gymnosperms, config$species$angiosperms)
  
  # Only keep the orders that are present in dt_res$order
  valid_orders <- all_orders[all_orders %in% dt_res$order]
  
  vebprint(valid_orders, verbose, "Valid Orders:")
  
  # Set the levels of the order factor
  dt_res[, order := factor(order, levels = valid_orders)]
  
  setorder(dt_res, order)
  
  vebprint(dt_res, verbose, "factor:")
  
  return(dt_res)
}

find_peaks <- function(data, column, threshold = 0.01, verbose = FALSE) {
  vebcat("Find Peaks Function:", veb = verbose)
  vebprint(data, verbose, "Input Data:")
  # Identify local maxima using diff
  peaks <- which(diff(sign(diff(data[[column]]))) < 0) + 1
  
  vebprint(peaks, verbose, "peaks:")
  
  # Filter based on threshold (difference from neighbors)
  if (!is.null(threshold)) {
    filtered_peaks <- peaks[data[[column]][peaks] - data[[column]][peaks - 1] >= threshold & data[[column]][peaks] - data[[column]][peaks + 1] >= threshold]
  } else {
    filtered_peaks <- peaks
  }
  
  vebprint(filtered_peaks, verbose, "Filtered Peaks:")
  
  out <- data[filtered_peaks, ]
  
  vebprint(out, text = "Out Data:", veb = verbose)
  
  # Return a data frame that only includes the peaks
  return(out)
}