get_sp_keys <- function(sp_names, out.dir, verbose = FALSE) {
  if (file.exists(paste0(out.dir, "/sp-w-keys.csv"))) {
    catn("File with keys found, loading file...")
    sp_w_keys <- fread(paste0(out.dir, "/sp-w-keys.csv"), sep = ",")
    return(sp_w_keys)
  }
  
  vebcat("Initiating GBIF Species keys download.", color = "funInit")
  
  create_dir_if(out.dir)
  
  if (!is.data.table(sp_names) && !is.vector(sp_names))  {
    vebcat("Error input is not a data frame or vector.", color = "fatalError")
    vebprint(str(sp_names), text = "Found class:")
    stop("Check class of input object.")
  }
  
  if (is.data.table(sp_names)) sp_names <- sp_names[[1]]
  
  l_names <- length(sp_names)
  
  catn("Collecting", highcat(l_names), "species keys...")

  taxon_info <- rgbif::name_backbone_checklist(sp_names)
  
  sp_w_keys <- data.table(
    scientificName = taxon_info$scientificName, 
    speciesKey = taxon_info$speciesKey
  )
  
  # Get NA keys and remove them
  na_sp_keys <- sp_w_keys[is.na(scientificName) | scientificName == "" | is.na(speciesKey) | speciesKey == ""]
  
  sp_w_keys <- sp_w_keys[!is.na(scientificName) & scientificName != "" & !is.na(speciesKey) & speciesKey != ""]
  
  # Write NA if any
  if (nrow(na_sp_keys) > 0) {
    catn("Removed", highcat(nrow(na_sp_keys)), "species with NA or blank keys.")
    na_out_file <- paste0(out.dir, "/na-sp-keys.csv")
    
    catn("Writing out na keys to:", colcat(na_out_file, color = "output"))
    
    fwrite(na_sp_keys, na_out_file, bom = TRUE)
  } else {
    vebcat("No blank keys, nor NAs found.", color = "proSuccess")
  }

  out_file <- paste0(out.dir, "/sp-w-keys.csv")

  catn("Writing out gbif_species to:", colcat(out_file, color = "output"))

  fwrite(sp_w_keys, out_file, row.names = F, bom = T)

  mdwrite(
    config$files$post_seq_md,
    text = paste0(
      "Number of NA species keys: **", nrow(na_sp_keys), "**  ",
      "Number of species keys for download: **", nrow(sp_w_keys), "**"
    )
  )

  vebcat("GBIF Species keys download completed successfully.", color = "funSuccess")

  return(sp_w_keys)
}
