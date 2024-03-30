filter_boreal <- function(spec.known, dfs, column, verbose = FALSE) {
  
  vebcat("Initiating boreal filter protocol.", color = "funInit")
  
  ##############
  # Initialize
  ##############
  
  arctic_present <- spec.known$present
  arctic_absent <- spec.known$absent
  glonaf_species <- dfs$glonaf_species
  
  boreal_dir <- "./outputs/filter/boreal"
  boreal_file <- paste0(boreal_dir, "/wfo-one-uniq.csv")
  boreal_log <- paste0(boreal_dir, "/logs")
  chunk_dir <- paste0(boreal_dir, "/chunk")
  chunk_name <- "species"
  output <- paste0(chunk_dir, "/", chunk_name)
  
  create_dir_if(c(boreal_dir, boreal_log, chunk_dir))
  
  
#########################
# Collect boreal species
#########################
  
  shapefiles <- c(
    boreal = "./resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp"
  )
  
  regions <- import_regions(
    shapefiles, 
    boreal_log, 
    verbose = T
  ) 
  
  boreal_region <- regions$region.boreal[regions$region.boreal$BIOME == 6, ]
  
  check_coords_orientation(boreal_region)
  
  boreal_wkt <- vect_to_wkt(
    boreal_region, 
    out.file = boreal_log,
    min.x = TRUE, 
    max.x = TRUE
  )
  
  check_coords_orientation(boreal_wkt)
  
  gbif_sp_list <- get_sp_list(
    taxon = "Tracheophyta", 
    region = boreal_wkt, 
    region.name = "boreal", 
    file.name = paste0(boreal_dir, "/boreal_sp_list"),
    download.key = "0037892-231120084113126", 
    download.doi = "https://doi.org/10.15468/dl.882wum"
  )
  
  # Run synonym Check on these names, then run the entire process again
  if (file.exists(boreal_file)) {
    catn("GBIF synonym check already conducted, reading file...")
    gbif_species <- fread(boreal_file)
  } else {
    filtered_gbif_list <- filter_gbif_list(gbif_sp_list = gbif_sp_list)
    
    # Got error when using "scientificName" in the WFO.match
    names(filtered_gbif_list) <- c("species", "rawName")
    
    gbif_sp_wfo <- check_syn_wfo(
      checklist = filtered_gbif_list, 
      column = "rawName", 
      folder = boreal_dir,
      max.cores = cores.max, 
      verbose = T, 
      counter = 100
    )
    
    gbif_sp_wfo_one <- check_syn_wfo_one(
      wfo_checklist = gbif_sp_wfo, 
      folder = boreal_dir
    )
    
    gbif_species <- gbif_sp_wfo_one$wfo_one_uniq
  }
  
  gbif_sp_list <- select_wfo_column(
    filepath = boreal_file,
    col.unique = "scientificName", 
    col.select = NULL,
    verbose = T
  )
  
  names(gbif_sp_list) <- "gbif_species"
  # Here we remove the glonaf species, because we do not need duplicates
  gbif_species <- dplyr::anti_join(gbif_sp_list$gbif_species, glonaf_species, by = column)
  
  
  ##################################################
  #           boreal arctic present absent         #
  ##################################################
  
  
  boreal_present <- write_filter_fun(
    file.out = paste0(boreal_dir, "/boreal-present-final.csv"),
    spec.in = gbif_species,
    spec.out = "Boreal present",
    fun = function() {
      # First merge to only get species from both dfs
      # This is because we do not want species not present in the arctic
      merged <- merge(gbif_species, arctic_present, by = column)
      
      return(merged)
    })
  
  # -------------------------------------------------------------------------- #
  
  boreal_absent <- write_filter_fun(
    file.out = paste0(boreal_dir, "/boreal-absent-final.csv"),
    spec.in = gbif_species,
    spec.out = "Boreal absent",
    fun = function() {
      # Remove species that are already present in the Arctic
      boreal_absent <-  dplyr::anti_join(gbif_species, arctic_present, by = column)
      
      # Update the absent list
      boreal_arc_absent = union_dfs(boreal_absent, arctic_absent, verbose = T)
      
      
      return(boreal_arc_absent)
    })
  
  vebcat("Boreal filter protocol completed successfully.", color = "funSuccess")
  
  return(list(
    spec = boreal_absent,
    dir = boreal_dir
  ))
}



