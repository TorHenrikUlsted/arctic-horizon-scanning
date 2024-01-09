get_sp_keys <- function(sp_names, out.dir, verbose = F) {
  cat(blue("Initiating GBIF Species keys download. \n"))
  
  if (file.exists(paste0(out.dir, "/sp-w-keys.csv"))) {
    cat("File with keys found, loading file... \n")
    sp_w_keys <- fread(paste0(out.dir, "/sp-w-keys.csv"), sep = ",")
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
      cat("found class:", class(sp_names), "\n")
      print(str(sp_names))
      stop("Check class of input object.")
    }
    
    sp_keys <- c()
    
    cat("Getting", cc$lightSteelBlue(length(sp_names)), "species keys. \n")
    
    cat("\ncollecting species keys \n\n")
    cat(sprintf("%6s | %11s | %8s\n", "key", "Total keys", "Percent"))  
    for (i in seq_along(sp_names)) {
      name <- sp_names[i]
      
      cat(cc$lightSteelBlue(sprintf("\r%6d | %11d | %8.2f %%", i, length(sp_names), round(i / length(sp_names) * 100, 2))))
      
      flush.console()
      
      taxon_info <- name_backbone(name)
      
      if (is.null(taxon_info) | !("usageKey" %in% names(taxon_info))) {
        sp_keys <- c(sp_keys, NA)
      } else if (is.na(taxon_info$usageKey)) {
        sp_keys <- c(sp_keys, NA)
      } else {
        sp_keys <- c(sp_keys, taxon_info$usageKey)
      }
    }
    
    sp_w_keys <- data.table(refinedScientificName = sp_names, usageKey = sp_keys)
    
    setorder(sp_w_keys, refinedScientificName)
    
    sp_w_keys <- set_df_utf8(sp_w_keys)
    
    cat("\n\nWriting out gbif_species to: \n", yellow("./outputs/setup/gbif/sp-w-keys.csv \n"))
    
    create_dir_if(out.dir)
    
    fwrite(sp_w_keys, paste0(out.dir, "/sp-w-keys.csv"), row.names = F, bom = T)
  }
  
  return(sp_w_keys)
}