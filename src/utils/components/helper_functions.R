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

#------------------------#
####      Spatial     ####
#------------------------#

is.spatVector <- function(x) {
  if (inherits(x, "SpatVector")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

is.spatRaster <- function(x) {
  if (inherits(x, "SpatRaster")) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

determine_data_nature <- function(x) {
  dt <- datatype(x)
  
  discrete_datatypes <- c("INT1U", "INT2S", "INT4S", "INT4U", "INT8S", "INT8U", "LOG1S")
  continuous_datatypes <- c("FLT4S", "FLT8S")
  
  if (is.spatVector(x)) {
    return("discrete")  # Always use "near" for SpatVectors
  }
  
  if (all(dt %in% continuous_datatypes)) {
    return("continuous")
  } else if (all(dt %in% discrete_datatypes)) {
    return("discrete")
  } else {
    warning("Mixed or unclear data types. Defaulting to 'near' method.")
    return("discrete")
  }
}

get_crs_config <- function(projection.name, vebose = FALSE) {
  
  try({
    projection <- crs(projection.name)
    return(projection)
  }, silent = TRUE)
  
  
  if (!is.character(projection.name)) {
    catn("Projection.name is not a string. Please chanche the input to a string.")
  }
  # Use config list to dynamically update the input_args if changing config values
  input_args <- config$projection$crs
  
  if (!projection.name %in% names(input_args)) {
    stop(paste("Invalid projection. Valid options are:", paste(names(input_args), collapse = ", ")))
  }
  
  if (inherits(crs(config$projection$crs$longlat, proj=T), "Error")) {
    cat("inehirts")
  }
  
  projection <- input_args[[projection.name]]
  
  catn("Using projection:", projection.name)
  
  return(projection)
}

extract_raster_to_dt <- function(raster, region = NULL, value = "value", cells = TRUE, xy = FALSE, verbose = FALSE) {
  vebcat("Extracting Raster and converting to data table.", veb = verbose)
  
  if (is.null(region)) {
    extract_by <- ext(raster)
  } else {
    extract_by <- region
  }
  
  rast_extr <- terra::extract(raster, extract_by, cells = cells, xy = xy)
  rast_dt <- as.data.frame(rast_extr)
  rast_dt <- as.data.table(rast_extr)
  
  if (is.null(region)) {
    names(rast_dt) <- c("cell", value)
  } else {
    if (xy) {
      names(rast_dt) <- c("ID", value, "cell", "longitude", "latitude")
    } else {
      names(rast_dt) <- c("ID", value, "cell")
    }
  }
  
  return(rast_dt)
}

convert_spatial_dt <- function(spatial, verbose = FALSE) {
  dt <- as.data.frame(spatial)
  dt <- as.data.table(dt)
  
  return(dt)
}

fix_antimeridian <- function(x, threshold = 0.00001, verbose) {
  catn("Fixing antimeridian issue.")
  
  ext_east <- terra::ext(ext(x)$xmin, 0, ext(x)$ymin, ext(x)$ymax)
  ext_west <- terra::ext(threshold, ext(x)$xmax, ext(x)$ymin, ext(x)$ymax)
  
  vect_east <- terra::crop(x, ext_east)
  vect_west <- terra::crop(x, ext_west)
  
  proj_east <- terra::project(vect_east, config$projection$crs$longlat)
  proj_west <- terra::project(vect_west, config$projection$crs$longlat)
  
  if (verbose) {
    plot(proj_east, main = "East side")
    plot(proj_west, main = "West side")
  }
  
  x_fixed <- rbind(proj_west, proj_east)
  
  if (verbose) plot(x_fixed, main = "Combined")
  
  return(x_fixed)
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

thin_occ_data <- function(dt, long = "decimalLongitude", lat = "decimalLatitude", projection = "+proj=longlat +datum=WGS84 +ellps=WGS84", res = 1000, seed = 123, verbose = FALSE) {
  setDT(dt)
  setalloccol(dt)
  
  # Add an ID to the input dt
  dt[, ID := .I]
  
  # get resolution in degrees
  resolution <- res / 111320
  
  # Make into point data
  sp_points <- vect(dt, geom = c(long, lat), crs = projection)
  vebprint(sp_points, verbose, "species points:")
  
  sp_ext <- ext(sp_points)
  
  # Create an empty raster with the same extent as the points
  r_points <- rast(sp_ext, crs = projection, res = resolution)
  values(r_points) <- 1:ncell(r_points)
  vebprint(r_points, verbose, "blank raster with same extent as points:")
  
  if (verbose) {
    plot(r_points)
    plot(sp_points, add = TRUE, col = "red")
  }
  
  # Extract the cell number of a point along with coordinates
  r_dt <- extract_raster_to_dt(r_points, sp_points)
  vebprint(r_dt, verbose, "Raster as data.table:")
  
  # Randomly keep unique points in each cell
  set.seed(seed)
  thinned_data <- r_dt[!is.na(cell), .SD[sample(.N, 1)], by = cell]
  
  # ignore NA cells because it could be points on the border, so do not thin them.
  vebprint(thinned_data, verbose, "Thinned data by choosing a random sample:")
  
  # Merge output with original data table
  thinned_dt <- dt[thinned_data, on = "ID", nomatch = 0]
  
  # Remove unnecessary columns
  thinned_dt[, c("ID", "cell", "value") := NULL]
  
  return(thinned_dt)
}

calc_lat_res <- function(lat_res, long_res, latitude = 0, unit.out = "km", verbose = FALSE) {
  lat_distance <- lat_res * 111.32
  
  long_deg_size <- 111.32 * cos(latitude * pi / 180)
  long_distance <- long_res * long_deg_size
  
  vebprint(lat_distance, text = "Latitude distance:")
  vebprint(long_distance, text = "Longitude distance:")
  
  # Find the highest resolution
  most_precise_res <- min(lat_distance, long_distance)
  
  vebprint(most_precise_res, text = "highest Resolution:")
  
  if (unit.out == "km") {
    return(most_precise_res)
  } else if (unit.out == "m") {
    return(floor(most_precise_res * 1000))
  }
}

calc_coord_uncertainty <- function(region, projection = "longlat", unit.out = "km", dir.out, verbose = FALSE) {
  out_file <- paste0(dir.out, "/coordinateUncertainty-", unit.out, ".txt")
  
  if (file.exists(out_file)) {
    max_res <- as.numeric(readLines(out_file))
  } else {
    catn("Calculating CoordinateUncertainty.")
    
    create_dir_if(dir.out)
    create_file_if(out_file)
    
    if (is.character(region)) {
      region <- rast(region)
    }
    
    if (terra::nlyr(region) > 1) {
      region <- terra::subset(region, 1)
    }
    
    region_ext <- terra::ext(region)
    
    vebprint(region_ext, text = "Region Extent:")
    
    region <- check_crs(region, projection = projection, projection.method = "bilinear")
    
    if (projection == "longlat") {
      res_lat <- terra::res(region)[2]
      res_long <- terra::res(region)[1]
      # Get latitude based on northern or southern hemisphere
      if (as.numeric(region_ext[4]) > 0) {
        lat <- as.numeric(region_ext[4]) # northern hemisphere
      } else {
        lat <- as.numeric(region_ext[3]) # Southern hemisphere
      }
      
      vebprint(lat, text = "Latitude:")
      
      max_res <- calc_lat_res(
        res_lat,
        res_long,
        lat,
        unit.out = unit.out,
        verbose = verbose
      )
    } else if (projection == "laea") {
      max_res <- floor(terra::res(region)[1])
      
      if (unit.out == "km") {
        max_res <- (max_res / 1000)
      }
    } else {
      stop("Error: only 'longlat' or 'laea' is available as projection parameters.")
    }
    
    catn("Writing file to:", colcat(out_file, color = "output"))
    
    writeLines(as.character(max_res), out_file)
  }
  
  catn("Lowest CoordinateUncertainty:", colcat(max_res, color = "indicator"))
  
  return(max_res)
}

load_sp_rast <- function(spec.filename) {
  sp_name <- basename(dirname(spec.filename))
  
  sp_rast <- terra::rast(spec.filename)
  names(sp_rast) <- sp_name
  
  return(sp_rast)
}

edit_crs <- function(crs.string, string.key, string.new, verbose = FALSE) {
  string.key <- toupper(string.key)
  
  # Split the keyword if it contains a number
  keyword_parts <- strsplit(string.key, "[[:digit:]]+", perl = TRUE)[[1]]
  keyword_num <- as.numeric(gsub("[^[:digit:]]", "", string.key))
  if (is.na(keyword_num)) keyword_num <- 1
  
  # Split the CRS string into parts
  crs_parts <- strsplit(crs.string, "\n", fixed = TRUE)[[1]]
  
  vebprint(crs_parts, verbose, "CRS parts:")
  
  # Find the parts that start with the keyword
  keyword_parts <- grep(paste0("\\b", keyword_parts, "\\b\\[\""), crs_parts)
  
  vebprint(keyword_parts, verbose, "Keyword parts:")
  
  # Check if the keyword exists in the CRS string
  if (length(keyword_parts) >= keyword_num) {
    # Replace the name following the keyword
    crs_parts[keyword_parts[keyword_num]] <- gsub("(?<=\\[\\\").*?(?=\\\")", string.new, crs_parts[keyword_parts[keyword_num]], perl = TRUE)
    
    # Combine the CRS parts back into a string
    new_crs <- paste(crs_parts, collapse = "\n")
  } else {
    catn("Keyword not found in the string.")
    # If the keyword doesn't exist, return the original CRS string
    new_crs <- crs.string
  }
  
  return(new_crs)
}

get_centroid_subregion <- function(region, region.sub = "subRegion", centroid.per.subregion = FALSE, inside = TRUE, verbose = FALSE) {
  uniq_subregions <- unique(region[[region.sub]])
  
  vebprint(uniq_subregions, verbose, "Unique Sub-Region(s):")
  
  if (centroid.per.subregion) {
    sub_region_centroids <- list()
  } else {
    sub_region_centroids <- vect()
  }
  
  for (i in 1:nrow(uniq_subregions)) {
    sub_region_name <- uniq_subregions[i, ]
    
    vebcat("Acquiring centroid for subregion", sub_region_name)
    
    vebprint(sub_region_name, verbose, "Sub-Region Name:")
    
    sub_region <- region[region[[region.sub]] == sub_region_name]
    
    vebprint(unique(sub_region[[region.sub]]), verbose, "Actual Subset Region:")
    
    vebprint(sub_region, verbose, "Sub Region:")
    
    all_centroids <- centroids(sub_region, inside = inside)
    
    vebprint(all_centroids, verbose, "All centroids:")
    
    n_centroids <- dim(all_centroids)[1]
    vebprint(n_centroids, verbose, "Dimensions:")
    
    if (n_centroids > 1) {
      all_x <- terra::crds(all_centroids)[, 1]
      all_y <- terra::crds(all_centroids)[, 2]
      
      mean_x <- mean(all_x)
      mean_y <- mean(all_y)
      
      euclidean_distances <- sqrt((all_x - mean_x)^2 + (all_y - mean_y)^2)
      
      centroid <- all_centroids[which.min(euclidean_distances), ]
    } else {
      centroid <- all_centroids
    }
    
    if (!centroid.per.subregion) {
      sub_region_centroids <- rbind(sub_region_centroids, centroid)
    } else {
      sub_region_centroids[[i]] <- centroid
      
      names(sub_region_centroids)[[i]] <- sub_region_name
    }
  }
  
  return(sub_region_centroids)
}
#------------------------#
####      ggplot      ####
#------------------------#

save_ggplot <- function(save.plot, save.name, save.width, save.height, save.dir, save.device = "jpeg", save.unit = "px", vis.title = FALSE, plot.show = FALSE, verbose = FALSE) {
  vebprint(save.plot, veb = plot.show)
  
  fig_out <- paste0(save.dir, "/", save.name, ".", save.device)
  
  create_dir_if(dirname(fig_out))
  
  catn("Saving plot to:", colcat(fig_out, color = "output"))
  ggsave(fig_out, device = save.device, unit = save.unit, width = save.width, height = save.height, plot = save.plot)
}

dynamic_guide_legend <- function(guide_config) {
  args <- list()
  
  for (param_name in names(guide_config)) {
    
    # Add the parameter to the args list
    args[[param_name]] <- guide_config[[param_name]]
  }
  
  # Create and return the guide_legend call
  do.call(guide_legend, args)
}

ggplot.filler <- function(gradient = "viridis-B", scale.type = "fill-c", limits = NULL, breaks = NULL, labels = NULL, begin = NULL, end = NULL, trans = NULL, guide = NULL, na.value = "transparent") {
  tryCatch(
    {
      # Syntax is: "gradient-option"
      split_str <- str_split(gradient, "-")[[1]]
      gradient <- split_str[[1]]
      option <- toupper(split_str[[2]])
      
      # Syntax is: "type-variable"
      split_str <- str_split(scale.type, "-")[[1]]
      scale_type <- split_str[[1]]
      scale_var <- tolower(split_str[[2]])
      
      # Process the guide
    if (!is.null(guide)) {
      if (is.list(guide)) {
        guide_obj <- dynamic_guide_legend(guide)
      } else if (is.character(guide) || is.logical(guide)) {
        guide_obj <- guide
      } else {
        warning("Invalid guide specification. Using default guide.")
        guide_obj <- "legend"
      }
    }
      
      args <- list(
        option = option,
        guide = if (!is.null(guide)) guide_obj else NULL,
        na.value = na.value
      )
      
      if (!is.null(labels)) {
        args$labels <- labels
      }
      
      if (!is.null(limits)) {
        args$limits <- limits
      }
      
      if (!is.null(breaks)) {
        args$breaks <- breaks
      }
      
      if (!is.null(begin)) {
        args$begin <- begin
      }
      
      if (!is.null(end)) {
        args$end <- end
      }
      
      if (!is.null(trans)) {
        args$trans <- trans
      }
      
      if (gradient == "viridis") {
        fun <- paste0("scale_", scale_type, "_viridis_", scale_var)
        return(do.call(fun, args))
      } else if (gradient == "whitebox") {
        fun <- paste0("scale_", scale_type, "_whitebox_", scale_var)
        args$palette <- args$option
        args$option <- NULL
        return(do.call(fun, args))
      }
    },
    error = function(e) {
      vebcat("Error when trying to use custom ggplot.filler function.", color = "fatalError")
      stop(e)
    }
  )
}

ggplot.theme <- function() {
  return(
    ggplot2::theme(
      text = element_text(family = config$ggplot$theme$text$family),
      
      plot.title = element_text(
        color = config$ggplot$theme$plot.title$color,
        vjust = config$ggplot$theme$plot.title$vjust,
        hjust = config$ggplot$theme$plot.title$hjust,
        size = config$ggplot$theme$plot.title$size,
        face = config$ggplot$theme$plot.title$face,
        margin = margin(b = config$ggplot$theme$plot.title$margin$b),
        lineheight = config$ggplot$theme$plot.title$lineheight
      ),
      
      plot.title.position = config$ggplot$theme$plot.title.position,
      plot.margin = do.call(margin, c(config$ggplot$theme$plot.margin[c("t", "r", "b", "l")], unit = config$ggplot$theme$plot.margin$unit)),
      
      axis.text = element_text(size = config$ggplot$theme$axis.text$size),
      axis.title.x = element_text(
        size = config$ggplot$theme$axis.title.x$size,
        hjust = config$ggplot$theme$axis.title.x$hjust
      ),
      
      axis.title.y = element_text(size = config$ggplot$theme$axis.title.y$size),
      legend.text = element_text(size = config$ggplot$theme$legend.text$size),
      legend.title = element_text(
        size = config$ggplot$theme$legend.title$size,
        hjust = config$ggplot$theme$legend.title$hjust
      ),
      
      legend.position = config$ggplot$theme$legend.position
    )
  )
}

#------------------------#
####      Objects     ####
#------------------------#

get_obj_name <- function(...) {
  sapply(as.list(match.call())[-1], deparse)
}

get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#------------------------#
####    Citations     ####
#------------------------#

format_bibtex <- function(bibtex_string) {
  # Remove leading/trailing whitespace and newline characters
  bibtex_string <- trimws(bibtex_string)
  bibtex_string <- gsub("\n", "", bibtex_string)
  
  # Split the string into individual fields
  fields <- strsplit(bibtex_string, ",(?=[^}]*(?:\\{|$))", perl = TRUE)[[1]]
  
  # Extract the entry type and key
  entry_type <- sub("^@(\\w+)\\{.*", "\\1", fields[1])
  
  # Process each field
  formatted_fields <- c(paste0("@", entry_type, "{,"))
  for (i in 2:length(fields)) {
    field <- fields[i]
    # Extract field name and value
    parts <- strsplit(field, "=")[[1]]
    field_name <- trimws(parts[1])
    field_value <- trimws(paste(parts[-1], collapse = "="))
    
    # Remove surrounding braces or quotes from field value
    field_value <- gsub("^[{\"](.*)[}\"]$", "\\1", field_value)
    
    # Format the field
    if (i == length(fields)) {
      # Last field: remove trailing } if present
      field_value <- sub("\\s*}\\s*$", "", field_value)
      formatted_fields <- c(formatted_fields, paste0("  ", field_name, " = {", field_value, "}"))
    } else {
      formatted_fields <- c(formatted_fields, paste0("  ", field_name, " = {", field_value, "},"))
    }
  }
  
  # Add the closing brace as a separate element
  formatted_fields <- c(formatted_fields, "}")
  
  # Return as a character vector
  return(formatted_fields)
}

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
  if (Sys.info()["sysname"] == "Windows") {
    total_rows <- as.numeric(system2("findstr", args = c("/R", "/N", "^", file.path), stdout = TRUE, stderr = NULL))
  } else { # for Unix-based systems like Linux and macOS
    total_rows <- as.numeric(system(paste("awk 'END {print NR}' ", file.path), intern = TRUE))
  }
  # -1 for header
  return(total_rows - 1)
}

system_calc_unique <- function(file.path, column, sep = "\t") {
  
  if (is.character(column)) {
    tmp <- fread(file.path, nrows = 0)
    column_number <- which(names(tmp) == column)
    if (length(column_number) == 0) stop("Column name not found")
  } else if (is.numeric(column)) {
    column_number <- as.integer(column)
  } else {
    stop("Column must be either a name (character) or a number")
  }
  
  if (Sys.info()["sysname"] == "Windows") {
    # Windows approach using R
    unique_values <- new.env(hash = TRUE)
    con <- file(file.path, "r")
    while (length(line <- readLines(con, n = 1)) > 0) {
      fields <- strsplit(line, sep, fixed = TRUE)[[1]]
      if (column_number <= length(fields)) {
        unique_values[[fields[column_number]]] <- TRUE
      }
    }
    close(con)
    unique_count <- length(unique_values)
  } else { # for Unix-based systems like Linux and macOS
    # Escape the separator for use in awk
    awk_sep <- gsub("\\", "\\\\", sep, fixed = TRUE)
    cmd <- sprintf("awk -F'%s' '{print $%d}' '%s' | sort -u | wc -l", awk_sep, column_number, file.path)
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
  
  # Escape the separator for use in awk
  awk_sep <- gsub("\\", "\\\\", sep, fixed = TRUE)
  
  if (Sys.info()["sysname"] == "Windows") {
    # Windows approach using PowerShell
    ps_command <- sprintf("Get-Content '%s' | ForEach-Object { $_.Split('%s')[%d] } | Group-Object | Measure-Object -Property Count, Name", 
                          file.path, sep, column_number - 1)
    result <- shell(sprintf("powershell -Command \"%s\"", ps_command), intern = TRUE)
    
    # Parse the result -1 for header
    total_rows <- as.numeric(strsplit(result[1], "\\s+")[[1]][2])
    unique_count <- as.numeric(strsplit(result[2], "\\s+")[[1]][2])
  } else {
    # Unix-based systems (Linux and macOS)
    cmd <- sprintf("awk -F'%s' '{print $%d}' '%s' | awk '{count++; unique[$0]++} END {print count, length(unique)}'", 
                   awk_sep, column_number, file.path)
    result <- system(cmd, intern = TRUE)
    
    # Parse the result
    values <- as.numeric(strsplit(result, " ")[[1]])
    total_rows <- values[1]
    unique_count <- values[2]
  }
  
  # -1 for header
  return(list(total_rows = total_rows - 1, unique_count = unique_count - 1))
}

model_to_md <- function(model) {
  output <- cli::cli_format({
    # Turn off CLI colors temporarily
    old <- options(cli.num_colors = 1)
    on.exit(options(old))
    
    # Execute the function
    output <- capture.output(summary(model), type = "output")
  })
  
  # Initialize empty string for markdown text
  md_text <- ""
  
  # Process the output line by line
  in_coefficients <- FALSE
  current_section <- NULL
  
  for (line in output) {
    # Remove quotes and leading/trailing whitespace
    line <- gsub('"', "", trimws(line))
    
    # Skip empty lines and divider lines
    if (grepl("^\\s*$", line) || grepl("^[-*]{10,}", line)) next
    
    if (grepl("link function:", line)) {
      # Start new section
      section_name <- sub(" link function:.*", "", line)
      md_text <- paste0(md_text, "\n## ", section_name, " Parameters\n\n")
      # Add coefficient table header
      md_text <- paste0(
        md_text,
        "| Parameter | Estimate | Std. Error | t value | p-value |\n",
        "|-----------|-----------|------------|----------|----------|\n"
      )
      in_coefficients <- TRUE
    }
    # Handle coefficient lines
    else if (in_coefficients && grepl("^\\(Intercept\\)|^median|^[a-zA-Z]", line)) {
      # Check if this is actually a statistics line
      if (grepl("^No\\. of observations|^Degrees of Freedom|^Global Deviance|^AIC|^SBC", line)) {
        in_coefficients <- FALSE
        if (!grepl("Model Statistics", md_text)) {
          md_text <- paste0(md_text, "\n## Model Statistics\n\n")
        }
        clean_line <- gsub("\\s+", " ", line)
        md_text <- paste0(md_text, "- ", clean_line, "\n")
        next
      }
      
      # Skip significance code lines
      if (grepl("^Signif\\.", line)) next
      
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      if (length(parts) >= 5) {
        param <- parts[1]
        est <- parts[2]
        se <- parts[3]
        t_val <- parts[4]
        
        # Handle p-value and stars
        p_val <- parts[5]
        if (length(parts) > 5) {
          if (grepl("\\*", p_val)) {
            p_val <- gsub("\\*", "\\\\*", p_val)
          } else {
            stars <- gsub("\\*", "\\\\*", parts[6])
            p_val <- paste0(p_val, " ", stars)
          }
        }
        
        md_text <- paste0(
          md_text,
          "| ", param, " | ",
          est, " | ",
          se, " | ",
          t_val, " | ",
          p_val, " |\n"
        )
      }
    } else if (grepl("^No\\. of observations|^Degrees of Freedom|^Global Deviance|^AIC|^SBC|^Residual Deg\\.|^at cycle:", line)) {
      if (!grepl("Model Statistics", md_text)) {
        md_text <- paste0(md_text, "\n## Model Statistics\n\n")
      }
      # Clean up the statistics line formatting
      clean_line <- gsub("\\s+", " ", line) # Replace multiple spaces with single space
      md_text <- paste0(md_text, "- ", clean_line, "\n")
    }
  }
  
  return(md_text)
}


mdwrite <- function(source, text = NULL, data = NULL, image = NULL, image.out = "./outputs/images/image", image.which = NULL, device = "jpeg", open = "a", veb = TRUE) {
  if (!veb) {
    return(invisible())
  }
  
  create_file_if(source, keep = TRUE)
  
  if (is.data.table(data) || is.data.frame(data)) {
    data <- kable(data, format = "markdown")
  }
  
  if (grepl(";", text)) {
    split_str <- str_split(text, ";")[[1]]
    h_num <- split_str[[1]]
    h_text <- split_str[[2]]
  }
  
  if (!is.null(image)) {
    if (grepl("\\.", basename(image.out))) {
      image_base <- sub("\\..*$", "", basename(image.out))
      image.out <- paste0(dirname(image.out), "/", image_base)
    }
    
    plot_title <- basename(image.out)
    image.out <- paste0(image.out, ".", device)
    
    catn("Writing image to:", colcat(image.out, color = "output"))
    
    do.call(device, list(filename = image.out, width = 500, height = 500))
    
    if (!is.null(image.which)) {
      plot(image, which = image.which, pch = "*", cex = 2, main = h_text)
    } else {
      plot(image, pch = "*", cex = 2, main = h_text)
    }
    
    dev.off()
  }
  
  try(con <- file(source, open = open))
  sink(con, type = "output")
  
  if (!is.null(text)) {
    if (grepl(";", text)) catn(paste0(strrep("#", h_num), " ", h_text))
    if (!grepl(";", text)) catn(text)
  }
  if (!is.null(data)) print(data)
  if (!is.null(image)) catn(paste0("![", h_text, "]", "(", "../images/", basename(image.out), ")"))
  catn("  ")
  
  sink(type = "output")
  close(con)
}

write_wrangled_md <- function(dt.list, name, column = "scientificName") {
  fl <- length(unique(dt.list$formatted, by = column)[[column]])
  
  counts <- lapply(dt.list, function(dt) {
    if (!is.null(dt)) length(unique(dt, by = column)[[column]])
  })
  
  sum <- sum(unlist(counts[names(counts) != "formatted"]), na.rm = TRUE)
  
  md_dt <- data.table(formatted = fl, lost = (fl - sum))
  
  for (key in names(dt.list)) {
    if (key != "formatted" && !is.null(dt.list[[key]])) {
      md_dt[[key]] <- counts[[key]]
    }
  }
  
  setcolorder(md_dt, c(setdiff(names(md_dt), "lost"), "lost"))
  
  mdwrite(
    config$files$post_seq_md,
    text = paste0("2;", name),
    data = md_dt
  )
}

create_derived_dataset <- function(occurrences.dir, verbose = FALSE) {
  sp_occ_out <- "./outputs/post-process/derived-data/datasetKey-count.csv"
  derived_data_zip_out <- "./outputs/post-process/derived-data/derived-dataset.zip"
  
  create_dir_if(dirname(sp_occ_out))
  
  if (file.exists(sp_occ_out)) {
    catn("DatasetKey with occurrence count already exists.")
    
    sp_occ_dt <- fread(sp_occ_out)
    
    catn(highcat(nrow(sp_occ_dt)), "datasetKeys added to csv file with corresponding total count.")
    
    catn(highcat(sum(sp_occ_dt$count)), "Total occurrences.")
    
    rm(sp_occ_dt)
    invisible(gc())
  } else {
    vebcat("Collecting datasetKeys and corresponding occurrence counts", color = "funInit")
    
    sp_occ_files <- list.files(occurrences.dir, full.names = TRUE)
    
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
    files <- list.files(occurrences.dir, full.names = FALSE)
    vebcat("Found", highcat(length(files)), "species in the directory.")
    rm(files)
    invisible(gc())
  } else {
    files <- list.files(occurrences.dir, full.names = TRUE)
    
    catn("Zipping", highcat(length(files)), "dervied species files.")
    
    zip(
      zipfile = derived_data_zip_out, 
      files = files,
      flags = "-j" # remove directories and only keep files
    )
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
