get_sp_keys <- function(sp_names, out.dir, verbose = F) {
  vebcat("Initiating GBIF Species keys download.", color = "funInit")
  
  if (file.exists(paste0(out.dir, "/sp-w-keys.csv"))) {
    catn("File with keys found, loading file...")
    sp_w_keys <- fread(paste0(out.dir, "/sp-w-keys.csv"), sep = ",")
  } else {
    if (!is.vector(sp_names)) {
      catn("Input failed class check, attempting to fix...")
      
      if(is.data.frame(sp_names)) {
        catn("Found a data frame, making into vector...")
        sp_names <- unlist(sp_names)
        if(is.vector(sp_names)) {
          vebcat("Successfully made data frame into vector.", color = "proSuccess")
        } else {
          stop("Failed to make data frame into vector.")
        }
      } 
      
    } else {
      vebcat("Error input is not a data frame or vector.", color = "fatalError")
      vebprint(str(sp_names), text = "Found class:")
      stop("Check class of input object.")
    }
    
    sp_keys <- c()
    
    catn("Getting", highcat(length(sp_names)), "species keys.")
    
    catn("\ncollecting species keys \n")
    cat(sprintf("%6s | %11s | %8s\n", "key", "Total keys", "Percent"))  
    for (i in seq_along(sp_names)) {
      name <- sp_names[i]
      
      cat(highcat(sprintf("\r%6d | %11d | %8.2f %%", i, length(sp_names), round(i / length(sp_names) * 100, 2))))
      
      flush.console()
      
      taxon_info <- name_backbone(name)
      
      if (is.null(taxon_info) | !("usageKey" %in% names(taxon_info))) {
        sp_keys <- c(sp_keys, NA)
      } else if (is.na(taxon_info$usageKey)) {
        sp_keys <- c(sp_keys, NA)
      } else {
        sp_keys <- c(sp_keys, taxon_info$usageKey)
      }
    }; catn()
    
    sp_w_keys <- data.table(species = sp_names, usageKey = sp_keys)
    
    setorder(sp_w_keys, species)
    
    # Remove and write out NAs
    keys_dt <- copy(sp_w_keys)
    na_keys_dt <- data.table()
    
    if (any(is.na(sp_w_keys$usageKey))) {
      catn("NA keys found, removing", highcat(length(which(is.na(sp_w_keys$usageKey)))), "species.")
      na_sp_keys <- keys_dt[is.na(keys_dt$usageKey), ]
      na_keys_dt <- rbind(na_keys_dt, na_sp_keys)
      
      sp_w_keys <- sp_w_keys[!is.na(sp_w_keys$usageKey)]
    } else if (any(sp_w_keys$usageKey) == "") {
      catn("Blank keys found, removing", highcat(length(which(sp_w_keys$usageKey == ""))))
      blnk_sp_keys <- keys_dt[keys_dt$usageKey == "", ]
      na_keys_dt <- rbind(na_keys_dt, blnk_sp_keys)
      
      sp_w_keys <- sp_w_keys[sp_w_keys$usageKey != ""]
    } else {
      vebcat("No blank keys, nor NAs found.", color = "proSuccess")
    }
    
    na_out_file <- paste0(out.dir, "/na-sp-keys.csv")
    
    catn("Writing out na keys to:", colcat(na_out_file, color = "output"))
    
    fwrite(na_keys_dt, na_out_file, bom = TRUE)
    
    sp_w_keys <- set_df_utf8(sp_w_keys)
    
    create_dir_if(out.dir)
    
    out_file <- paste0(out.dir, "/sp-w-keys.csv")
    
    catn("Writing out gbif_species to:", colcat(out_file, color = "output"))
    
    fwrite(sp_w_keys, out_file, row.names = F, bom = T)
    
    mdwrite(
      post_seq_nums,
      heading = paste0(
        "Number of NA species keys: **", nrows(na_keys_dt), "**  ",
        "Number of species keys for download: **", nrow(sp_w_keys), "**"
      )
    )
    
  }
  
  vebcat("GBIF Species keys download completed successfully.", color = "funSuccess")
  
  return(sp_w_keys)
}