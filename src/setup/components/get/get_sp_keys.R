get_sp_keys <- function(sp_names, out.dir, verbose = F) {
  cat(blue("Initiating GBIF Species keys download. \n"))
  
  if (file.exists(paste0(out.dir, "/sp_w_keys.csv"))) {
    cat("File with keys found, loading file... \n")
    sp_w_keys <- fread(paste0(out.dir, "/sp_w_keys.csv"), sep = ",")
  } else {
    if (!is.vector(sp_names)) {
      cat(cc$lightCoral("Input failed class check, attempting to fix... \n"))
      
      if(is.data.frame(sp_names)) {
        cat("Found a data frame, making into vector... \n")
        sp_names <- unlist(sp_names)
        if(is.vector(sp_names)) cat(cc$lightGreen("Successfully made data frame into vector. \n")) else stop("Failed to make data frame into vector.")
      }
      
    } else {
      cat(red("Error input is not a data frame or vector. \n"))
      stop("Check input")
    }
    
    sp_keys <- c()
    
    cat("Getting", cc$lightSteelBlue(length(sp_names)), "species keys. \n")
    
    cat("\ncollecting species keys \n\n")
    cat(sprintf("%6s | %11s | %8s\n", "key", "Total keys", "Percent"))  
    for (i in seq_along(sp_names)) {
      name <- sp_names[i]
      
      cat(sprintf("\r%6d | %11d | %8.2f %%", cc$lightSteelBlue(i), cc$lightSteelBlue(length(sp_names)), cc$lightSteelBlue(round(i / length(sp_names) * 100, 2))))
      
      flush.console()
      
      taxon_info <- name_backbone(name)
      
      if (is.na(taxon_info$speciesKey)) {
        sp_keys <- c(sp_keys, NA)
      } else {
        sp_keys <- c(sp_keys, taxon_info$speciesKey)
      }
    }
    
    sp_w_keys <- data.table(scientificName = sp_names, speciesKey = sp_keys)
    
    sp_w_keys <- set_df_utf8(sp_w_keys)
    
    cat("\nWriting out gbif_species to: \n", yellow("./outputs/setup/gbif/sp_w_keys.csv \n"))
    
    create_dir_if(out.dir)
    
    fwrite(sp_w_keys, paste0(out.dir, "/sp_w_keys.csv"), row.names = F, bom = T)
  }
  
  return(sp_w_keys)
}