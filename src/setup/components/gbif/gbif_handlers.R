gbif.match <- function(checklist) {
  sp_w_keys <- as.data.table(gbif_retry(checklist$verbatimName.GBIF, "name_backbone_checklist"))
  # Handle cases where some columns are not returned by GBIF
  cols <- c("verbatim_name", "species", "scientificName", "rank", "speciesKey", "usageKey", "acceptedUsageKey", "status", "synonym")
  
  # Ensure all expected columns exist, filling with NA if missing
  for (col in cols) {
    if (!(col %in% names(sp_w_keys))) {
      sp_w_keys[, (col) := NA]
    }
  }
  
  # Subset columns
  sp_w_keys <- sp_w_keys[, ..cols]
  
  old_names <- cols[-(2:3)]
  new_names <- c("verbatimName.GBIF", "rank.GBIF", "speciesKey.GBIF", "usageKey.GBIF", "acceptedUsageKey.GBIF", "status.GBIF", "synonym.GBIF")
  
  data.table::setnames(sp_w_keys, old_names, new_names) 
  
  return(sp_w_keys)
}

gbif.accepted <- function(match.result) {
  for (i in 1:nrow(match.result)) {
    cat("\rChecking scientificName", highcat(i), "/", nrow(match.result))
    
    spec <- match.result[i]
    
    if (spec$synonym.GBIF == FALSE) next
    
    verbatimName <- spec$verbatimName.GBIF
    search_key <- spec$acceptedUsageKey.GBIF
    
    dt <- as.data.table(gbif_retry(search_key, "name_usage")$data)
    
    # vebprint(match.result[i], text = "match.res:")
    # vebprint(dt, text = "dt:")
    # 
    # catn(verbatimName)
    # catn(dt$species)
    # catn(dt$scientificName)
    # catn(dt$rank)
    # catn(dt$speciesKey)
    # catn(dt$key)
    # catn(dt$taxonomicStatus)
    
    match.result[i] <- data.table(
      verbatimName.GBIF = verbatimName,
      species = if (!is.na(match.result[i]$species)) dt$species else NA,
      scientificName = dt$scientificName,
      rank.GBIF = dt$rank,
      speciesKey.GBIF = if (!is.na(match.result[i]$speciesKey.GBIF)) dt$speciesKey else NA,
      usageKey.GBIF = dt$key,
      acceptedUsageKey.GBIF = NA,
      status.GBIF = dt$taxonomicStatus,
      synonym.GBIF = FALSE
    )
  };catn()
  
  return(match.result)
}

gbif_standardize <- function(dt, out.file, verbose = FALSE) {
  if (nrow(dt) == 0) {
    catn(highcat(nrow(dt)), "scientificNames found")
    return(dt)
  }
  
  if (file.exists(out.file)) return(fread(out.file))
  
  vebcat("Identifying accepted GBIF usageKeys", color = "indicator")
  
  orig_n <- nrow(dt)
  
  if (grepl("wfo", out.file)) {
    data.table::setnames(dt, c("scientificName", "scientificNameAuthorship"), c("scientificName.WFO", "scientificNameAuthorship.WFO"))
    dt[, `:=` (
      verbatimName.GBIF = fifelse(
        !is.na(scientificNameAuthorship.WFO) & scientificNameAuthorship.WFO != "", 
        paste(scientificName.WFO, scientificNameAuthorship.WFO), 
        scientificName.WFO
      )
    )]
  } else {
    dt[, `:=` (
      verbatimName.GBIF = scientificName, # Keep the original
      scientificName = NULL # Remove the column
    )]
  }
  
  dt <- unique(dt, by="verbatimName.GBIF")
  init_dups <- orig_n - nrow(dt)
  
  gbif_matched <- gbif.match(dt)
  
  gbif_accepted <- gbif.accepted(gbif_matched)
  
  gbif_accepted <- dt[gbif_accepted, on = "verbatimName.GBIF", nomatch = 0L]
  
  # Remove anything NA or names above species taxon level
  gbif_na <- gbif_accepted[is.na(scientificName) | scientificName == "" | is.na(speciesKey.GBIF) | speciesKey.GBIF == ""]
  
  na_n <- nrow(gbif_na)
  
  if (!grepl("mismatch", out.file)) {
    if (na_n > 0) {
      catn("Found", highcat(na_n), "NA speciesKeys or scientificNames")
      na_out_name <- basename(gsub(".csv", "", out.file))
      na_out_file <- file.path(dirname(out.file), paste0(na_out_name, "-na", ".csv"))
      
      catn("Writing out na names or speciesKeys to:", colcat(na_out_file, color = "output"))
      fwrite(gbif_na, na_out_file, bom = TRUE)
    }
    
    res_dt <- gbif_accepted[!is.na(scientificName) & scientificName != "" & !is.na(speciesKey.GBIF) & speciesKey.GBIF != ""]
    
    after_na_removal_n <- nrow(res_dt)
    
    res_dt <- unique(res_dt, by = "scientificName")
    unique_n <- nrow(res_dt)
  } else {
    res_dt <- gbif_accepted
    after_na_removal_n <- 0
    unique_n <- 0
  }
  
  new_n <- nrow(res_dt)
  
  dups_after <- after_na_removal_n - unique_n
  tot_dups <- init_dups + dups_after
  
  catn(highcat(tot_dups), "duplicate scientific names removed")
  
  if (grepl("wfo", out.file)) {
    changed <- length(which(paste(res_dt$verbatimName.GBIF, res_dt$verbatimNameAuthorship.GBIF) != res_dt$scientificName))
  } else {
    changed <- length(which(res_dt$verbatimName.GBIF != res_dt$scientificName))
  }
  
  md_dt <- data.table(
    input = orig_n,
    output = new_n,
    changes = changed,
    duplicate = tot_dups,
    na = na_n
  )
  
  mdwrite(
    config$files$post_seq_md,
    text = paste0("3;GBIF accepted name conversion for ", basename(dirname(out.file))),
    data <- md_dt
  )
  
  fwrite(res_dt, out.file, bom = TRUE)
  
  return(res_dt)
}

filter_keys <- function(keys.dt, out.dir, verbose = FALSE) {
  name <- basename(dirname(out.dir))
  
  data.table::setnames(keys.dt, "usageKey.GBIF", "usageKey")
  # Filter out usageKeys that will be included in speciesKey
  species_keys <- keys.dt[tolower(rank.GBIF) == "species", unique(speciesKey.GBIF)]
  # Filter out lower-level taxa if their species is already present
  removed_keys <- keys.dt[speciesKey.GBIF %in% species_keys & tolower(rank.GBIF) != "species"]
  
  if (nrow(removed_keys) > 0) {
    catn("Found", highcat(nrow(removed_keys)), "taxons below species where the species is already included")
    fwrite(removed_keys, paste0(out.dir, "/dup_sp_keys.csv"))
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0("**", nrow(removed_keys), "** taxons below an already included species were removed for ", name, " species")
    )
  }
  
  keys_dt <- keys.dt[!(speciesKey.GBIF %in% species_keys & tolower(rank.GBIF) != "species")]
  
  na_keys <- keys.dt[(is.na(usageKey))]
  
  if (nrow(na_keys) > 0) {
    catn("Found", highcat(nrow(na_keys)), "NA usageKeys when filtering keys")
    fwrite(na_keys, paste0(out.dir, "/na_sp_keys.csv"))
    
    mdwrite(
      config$files$post_seq_md,
      text = paste0("**", nrow(na_keys), "** NA usageKeys were removed for ", name, " species")
    )
  }
  
  na_keys <- keys.dt[!(is.na(usageKey))]
  
  fwrite(keys_dt, paste0(out.dir, "/filtered-sp-keys.csv"))
  
  keys_dt <- keys_dt[, .(usageKey)]
  
  rm(keys.dt, name, species_keys, removed_keys)
  invisible(gc())
  
  return(keys_dt)
}
