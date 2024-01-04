source_all("./src/filter/components")

filter = function(verbose) {
  cat(blue("Initiating filtering sequence. \n"))
  
  
  ##################################################
  #                     setup                      #
  ##################################################
  
  
  check_cpu_speed(
    df.path = "./resources/data-raw/cpu-test-species.csv", 
    max.cores = 1,
    verbose = T,
    counter = 1
  )
  
  filter_timer = start_timer("filter_timer")
  
  setup_raw_data(
    column = "rawName", 
    test = NULL, 
    max.cores = 20,
    verbose = T,
    counter = 1
    )
  

  ##################################################
  #                     load dfs                   #
  ##################################################
  
  if (verbose) cat("Loading dfs. \n")
  
  dfs <- select_wfo_column(read.dir = "./resources/synonym-checked", column = "scientificName", verbose = T)
  
  # Combine no-matches
  man_checked_nomatch <- fread("./outputs/setup/wrangle/wfo-nomatch-edited.csv")
  
  # Remove NAs and blank values in the newName df
  man_formatted <- man_checked_nomatch %>% 
    dplyr::filter(!is.na(newName) & newName != "") %>% 
    dplyr::select(newName, listOrigin)
  
  # Set them to scientificNames to be the same as the others
  names(man_formatted)[1] <- "scientificName"
  
  for(i in 1:nrow(man_formatted)){
    # Get the scientificName and listOrigin from the current row
    scientificName <- man_formatted[i, "scientificName"]
    listOrigin <- as.character(man_formatted[i, "listOrigin"])
    
    cat("Adding", cc$lightSteelBlue(scientificName), "to", cc$lightSteelBlue(listOrigin), "\n")
    
    # Check if scientificName or listOrigin is missing
    if (is.na(scientificName) | is.na(listOrigin)) {
      cat("Missing value in row ", i, "\n")
      next
    }
    
    # Check if listOrigin exists in dfs
    if (!listOrigin %in% names(dfs)) {
      cat("Invalid listOrigin in row ", i, ": ", listOrigin, "\n")
      next
    }
    
    # Append the scientificName to the corresponding data table in dfs
    dfs[[listOrigin]] <- rbindlist(list(dfs[[listOrigin]], data.table(scientificName)))
  }
  
  cat(cc$lightGreen("the manually formatted synonym checks have been successfully added to correct data frames. \n"))
  
  ##################################################
  #              combine aba ambio                 #
  ##################################################
  
  aba_present = dfs$aba_present
  ambio_present  = dfs$ambio_present
  
  arctic_present <- union_dfs(aba_present, ambio_present, verbose = T)
  
  aba_absent = dfs$aba_absent
  ambio_absent  = dfs$ambio_absent
  
  arctic_absent <- union_dfs(aba_absent, ambio_absent, verbose = T)
  
  ##################################################
  #         glonaf arctic present absent           #
  ##################################################
  
  glonaf_species <- dfs$glonaf_species
  
  # First merge to only get species from both dfs
  glonaf_present <- merge(glonaf_species, arctic_present, by = "scientificName")
  
  cat(cc$lightSteelBlue(nrow(glonaf_present)), "GloNAF species already exist in the Arcitc. \n")
  
  glonaf_absent <-  anti_join(glonaf_species, arctic_present, by = "scientificName")
  
  # Update the absent list
  arctic_absent = union_dfs(glonaf_absent, arctic_absent, verbose = T)
  
  create_dir_if("./outputs/filter/arctic")
  fwrite(arctic_absent, "./outputs/filter/arctic/arctic-absent-final.csv", bom = T)
  
  ##################################################
  #               Collect gbif species             #
  ##################################################
  
  shapefiles <- c(
    boreal = "./resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp"
  )
  
  regions <- import_regions(
    shapefiles, 
    "./outputs/filter/logs", 
    verbose = T
    ) 
  
  boreal_region <- regions$boreal[regions$boreal$BIOME == 6, ]
  
  check_coords_orientation(boreal_region)
  
  boreal_wkt <- vect_to_wkt(boreal_region, min.x = T, max.x = T)
  
  check_coords_orientation(boreal_wkt)
  
  gbif_sp_list <- get_sp_list(
    taxon = "Tracheophyta", 
    region = boreal_wkt, 
    region.name = "boreal", 
    file.name = "./outputs/filter/gbif-acquisition/gbif-sp-list", 
    download.key = "0037892-231120084113126", 
    download.doi = "https://doi.org/10.15468/dl.882wum"
    )
  
  cat("Filtering the GBIF list. \n")
  
  gbif_filtered <- gbif_sp_list %>% 
    dplyr::filter(taxonRank %in% c("SPECIES", "SUBSPECIES", "VARIETY", "FORM", "UNRANKED"))
  
  gbif_sp_scientificName <- data.frame(scientificName = gbif_filtered$acceptedScientificName)
  
  if (any(is.na(gbif_sp_scientificName))) {
    cat(cc$lightCoral("Some species are NA, removing... \n"))
    gbif_sp_scientificName <- gbif_sp_scientificName[!is.na(rowSums(gbif_sp_scientificName)), ]
    
    if (any(is.na(gbif_sp_scientificName))) cat(cc$lightCoral("Failed to remove NA values. \n")) else cat(cc$lightGreen("NAs removed successfully. \n"))
  } else {
    cat(cc$lightGreen("No NAs found. \n"))
  }
  
  if (any(gbif_sp_scientificName == "")) {
    cat(cc$lightCoral("Some species are blank, removing... \n"))
    gbif_sp_scientificName <- gbif_sp_scientificName[rowSums(gbif_sp_scientificName != "") > 0, ]
    if (any(is.na(gbif_sp_scientificName))) cat(cc$lightCoral("Failed to remove NA values. \n")) else cat(cc$lightGreen("NAs removed successfully. \n"))
  } else {
    cat(cc$lightGreen("No blank valuess found. \n"))
  }
  
  cat("Removing duplicates. \n")
  
  # Remove species that are already present in the Arctic
  gbif_sp_scientificName <-  anti_join(gbif_sp_scientificName, arctic_present, by = "scientificName")
  
  # Remove species already a part of the glonaf_absent species
  gbif_sp_scientificName <-  anti_join(gbif_sp_scientificName, arctic_absent, by = "scientificName")
  
  names(gbif_sp_scientificName) <- "rawName"
  
  # Run synonym Check on these names, then run the entire process again
  if (file.exists("./outputs/filter/gbif-acquisition/wfo-one-uniq.csv")) {
    cat("GBIF synonym check already conducted, reading file... \n")
    gbif_species <- fread("./outputs/filter/gbif-acquisition/wfo-one-uniq.csv")
  } else {
    gbif_sp_wfo <- check_syn_wfo(
      checklist = gbif_sp_scientificName, 
      column = "rawName", 
      folder = "./outputs/filter/gbif-acquisition", 
      max.cores = 20, 
      verbose = T, 
      counter = 100
    )
    
    gbif_sp_wfo_one <- check_syn_wfo_one(
      wfo_checklist = gbif_sp_wfo, 
      folder = "./outputs/filter/gbif-acquisition"
    )
    
    gbif_species <- gbif_sp_wfo_one$wfo_one_uniq
    
  }
  
  gbif_species <- gbif_species %>% 
    select(scientificName)
  
  
  if (any(duplicated(gbif_species))) {
    cat(cc$lightSteelBlue(length(which(duplicated(gbif_species)))), cc$lightCoral("Species are duplicated. Removing... \n"))
    
    gbif_species <- unique(gbif_species)
    if (any(duplicated(gbif_species))) cat(cc$lightCoral("Failed to remove duplicates. \n")) else cat(cc$lightGreen("Successfully removed duplicated species. \n"))
  } else {
    cat(cc$lightGreen("All species are unique. \n"))
  }
  
  ##################################################
  #           gbif arctic present absent           #
  ##################################################
  
  
  cat("Separating", cc$lightSteelBlue(nrow(gbif_species)), "gbif species into arctic present and absent. \n")
  
  # First merge to only get species from both dfs
  gbif_present <- merge(gbif_species, arctic_present, by = "scientificName")
  
  cat("GBIF has",  cc$lightSteelBlue(nrow(gbif_present)), "species already present in the Arctic. \n")
  
  # Remove species that are already present in the Arctic
  gbif_absent <-  anti_join(gbif_species, arctic_present, by = "scientificName")
  
  # Remove species already a part of the arctic_absent species
  gbif_arc_absent <-  anti_join(gbif_absent, arctic_absent, by = "scientificName")
  
  cat("Removed", cc$lightSteelBlue(nrow(gbif_absent) - nrow(gbif_arc_absent)), "GBIF species already existing in the Arctic absent list. \n")
  
  cat("There are", cc$lightSteelBlue(nrow(gbif_arc_absent)), "Unique GBIF species not in the Arctic absent list. \n")
  
  # Update the absent list
  #arctic_absent = union_dfs(gbif_glo_absent, arctic_absent, verbose = T)
  
  
  ##################################################
  #        Occurrence arctic absent download       #
  ##################################################
  
  arctic_absent_keys <- get_sp_keys(
    sp_names = arctic_absent,
    out.dir = "./outputs/filter/arctic",
    verbose = T
    )
  
  arctic_absent_occ <- get_occ_data(
    species_w_keys = arctic_absent_keys,
    file.name = "outputs/filter/arctic/arctic-absent-occ", 
    region = NULL, 
    download.key = "0042685-231120084113126",
    download.doi = "https://doi.org/10.15468/dl.bwafdv"
  )
  
  
  ##################################################
  #        Occurrence gbif absent download         #
  ##################################################
  
  # gbif_absent_keys <- get_sp_keys(
  #   sp_names = gbif_absent,
  #   out.dir = "./outputs/filter/gbif-occ",
  #   verbose = T
  #   )
  # 
  # gbif_absent_occ <- get_occ_data(
  #   species_w_keys = gbif_absent_keys,
  #   file.name = "outputs/occurrence/gbif/gbif-absent-occ", 
  #   region = NULL, 
  #   download.key = NULL, 
  #   download.doi = NULL
  # )
  
  end_timer(filter_timer)
}
