union_dfs <- function(df1, df2, verbose = F) {
  
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

fix_nomatches <- function(dfs, nomtach.edited, column) {
  # Combine no-matches
  edited_nomatch <- fread(nomatch.edited)
  
  # Remove NAs and blank values in the newName df
  formatted <- edited_nomatch %>% 
    dplyr::filter(!is.na(all_of(column)) & all_of(column) != "") %>% 
    dplyr::select(scientificNameAuthorship, all_of(column), listOrigin)
  
  vebprint(formatted, text = "Manually formatted data:")
  
  for(i in 1:nrow(formatted)){
    
    # Get the scientificName and listOrigin from the current row
    scientificNameAuthorship <- formatted[i, "scientificNameAuthorship"]
    scientificName <- formatted[i, ..column]
    listOrigin <- as.character(formatted[i, "listOrigin"])
    
    catn("Appending", highcat(as.character(scientificName)), "to", highcat(listOrigin))
    
    # Check if scientificName or listOrigin is missing
    if (is.na(scientificName) | is.na(listOrigin)) {
      catn("Missing value in row ", i)
      next
    }
    
    # Check if listOrigin exists in dfs
    if (!listOrigin %in% names(dfs)) {
      catn("Invalid listOrigin in row ", i, ": ", listOrigin)
      next
    }
    
    # Create a new data table with the same columns as dfs[[listOrigin]]
    new_row <- data.table(scientificName)
    setnames(new_row, names(dfs[[listOrigin]]))
    new_row[1, ] <- c(scientificName)
    
    # Append the scientificName to the corresponding data table in dfs
    dfs[[listOrigin]] <- rbindlist(list(dfs[[listOrigin]], new_row))
  }
  
  vebcat("the manually formatted synonym checks have been successfully added to correct data frames.", veb = verbose, color = "proSuccess")
  
  return(dfs)
}

write_filter_fun <- function(file.out, spec.in, spec.out, fun) {
  create_dir_if(dirname(file.out))
  
  catn(highcat(nrow(spec.in)), "Species input.")
  
  result <- fun
  
  catn(highcat(nrow(result)), "species output")
  catn("In/out difference:", highcat(nrow(spec.in) - nrow(result)))
  
  catn("Writing", spec.out, "to:", colcat(file.out, color = "output"))
  
  fwrite(result, file.out, bom = T)
  
  return(result)
}

# For file.out, do not use extension
get_occurrence <- function(spec, file.out, region, download.key, doi.key, verbose) {
  
  gbif_absent_keys <- get_sp_keys(
    sp_names = spec,
    out.dir = dirname(file.out),
    verbose = verbose
  )
  
  gbif_absent_occ <- get_occ_data(
    species_w_keys = gbif_absent_keys,
    file.name = file.out,
    region = region,
    download.key = download.key,
    download.doi = download.doi
  )
  
  return(occ_data)
}

# spec.occ can be either filepath or data frame
chunk_protocol <- function(
    spec.occ,
    chunk.name = "species",
    chunk.col,
    chunk.dir,
    chunk.size = 1e6,
    cores.max = 1,
    iterations = NULL,
    verbose
  ) {
  
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
      df = sp_occ,
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
    chunk.col = "combined"
  }
  
  clean_chunks(
    chunk.name = chunk.name,
    chunk.column = chunk.col, 
    chunk.dir = chunk.dir, 
    sp_w_keys = sp_w_keys_out,
    verbose = verbose
  )
}


