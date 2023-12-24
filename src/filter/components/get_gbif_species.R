source("./src/filter/components/get_gbif_keys_region.R")
source("./src/filter/components/get_gbif_names_keys.R")

get_gbif_species <- function(region, region.name, taxon, cores.max, verbose) {

  gbif_counts <- get_gbif_keys_region(region, region.name, taxon, verbose) 

  gbif_names <- get_gbif_names_keys(region, gbif_counts, cores.max, verbose)
  
  print(head(gbif_names, 10))
  
  scientific_names <- unlist(gbif_sp_wfo_one$scientificName)
  
  gbif_species_check <- rgbif_simliarity_check(scientific_names, region)
  
  print(class(gbif_species_check))
  print(str(gbif_species_check))
  print(head(gbif_species_check, 10))
  
  
  ## Add to top
  if (file.exists(check_savefile)) {
    cat("Occurrence data found, Loading... \n")
    
    gbif_species <- fread(check_savefile, sep = "\t")
    
  } else {
    
  }

  return(sp_occ_cropped)
}
