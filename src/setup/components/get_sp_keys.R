get_sp_keys <- function(sp_names) {
  cat(blue("Initiating GBIF Species keys download. \n"))
  
  sp_keys <- c()
  
  for (name in sp_names) {
    taxon_info <- name_backbone(name)
    
    sp_keys <- c(sp_keys, taxon_info$speciesKey)
  }
  
  sp_w_keys <- data.table(scientificName = sp_names, speciesKey = sp_keys)
  
  sp_w_keys <- set_df_utf8(sp_w_keys)
  
  cat(cyan("Writing out gbif_species to: \n", "./outputs/setup/gbif/sp_w_keys.csv \n"))
  
  if (!dir.exists(paste0("./outputs/setup/gbif/"))) dir.create(paste0("./outputs/setup/gbif/"), recursive = T)
  
  fwrite(sp_w_keys, "./outputs/setup/gbif/sp_w_keys.csv", row.names = F, bom = T)
  
  return(sp_w_keys)
}