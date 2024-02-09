source_all("./src/filter/components")

filter_process <- function(test = NULL, cores.max = 1, verbose) {
  on.exit(closeAllConnections())
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
  on.exit(end_timer(filter_timer))
  
  setup_raw_data(
    column = "rawName",
    test = test, 
    max.cores = cores.max,
    verbose = T,
    counter = 1
  )
  
  ##################################################
  #                     load dfs                   #
  ##################################################
  
  if (verbose) cat("Loading dfs. \n")
  
  dfs <- select_wfo_column(
    filepath = "./resources/synonym-checked", 
    col.unique = c("genus", "specificEpithet", "infraspecificEpithet"), 
    col.select = "scientificNameAuthorship",
    verbose = T
    )

  if (is.null(test)) {
    
    # Combine no-matches
    man_checked_nomatch <- fread("./outputs/setup/wrangle/wfo-nomatch-edited.csv")
    
    # Remove NAs and blank values in the newName df
    man_formatted <- man_checked_nomatch %>% 
      dplyr::filter(!is.na(refinedScientificName) & refinedScientificName != "") %>% 
      dplyr::select(refinedAuthor, refinedScientificName, listOrigin)
    
    print(man_formatted)
    
    for(i in 1:nrow(man_formatted)){
      # Get the scientificName and listOrigin from the current row
      scientificNameAuthorship <- man_formatted[i, "refinedAuthor"]
      refinedScientificName <- man_formatted[i, "refinedScientificName"]
      listOrigin <- as.character(man_formatted[i, "listOrigin"])
      
      cat("Appending", cc$lightSteelBlue(refinedScientificName), "to", cc$lightSteelBlue(listOrigin), "\n")
      
      # Check if scientificName or listOrigin is missing
      if (is.na(refinedScientificName) | is.na(listOrigin)) {
        cat("Missing value in row ", i, "\n")
        next
      }
      
      # Check if listOrigin exists in dfs
      if (!listOrigin %in% names(dfs)) {
        cat("Invalid listOrigin in row ", i, ": ", listOrigin, "\n")
        next
      }
      
      # Create a new data table with the same columns as dfs[[listOrigin]]
      new_row <- data.table(scientificNameAuthorship, refinedScientificName)
      setnames(new_row, names(dfs[[listOrigin]]))
      new_row[1, ] <- c(scientificNameAuthorship, refinedScientificName)
      
      #cat("Using column", cc$lightSteelBlue(names(dfs[[listOrigin]])), "for table 1 and", cc$lightSteelBlue(names(data.table(new_row))), "for table 2. \n")
      # Append the scientificName to the corresponding data table in dfs
      dfs[[listOrigin]] <- rbindlist(list(dfs[[listOrigin]], new_row))
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
    glonaf_present <- merge(glonaf_species, arctic_present, by = "refinedScientificName")
    
    cat(cc$lightSteelBlue(nrow(glonaf_present)), "GloNAF species already exist in the Arcitc. \n")
    
    glonaf_absent <-  anti_join(glonaf_species, arctic_present, by = "refinedScientificName")
    
    # Update the absent list
    arctic_absent = union_dfs(glonaf_absent, arctic_absent, verbose = T)
    
    create_dir_if("./outputs/filter/arctic")
    
    arctic_absent <- arctic_absent[, c("refinedScientificName", "scientificNameAuthorship", setdiff(names(arctic_absent), c("refinedScientificName", "scientificNameAuthorship"))), with = FALSE]
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
    
    gbif_sp_scientificName <- data.table(species = gbif_filtered$species, scientificName = gbif_filtered$scientificName)
    
    if (any(is.na(gbif_sp_scientificName))) {
      cat(cc$lightCoral("Some species are NA, removing... \n"))
      gbif_sp_scientificName <- data.table(gbif_sp_scientificName[!is.na(rowSums(gbif_sp_scientificName)), ])
      
      if (any(is.na(gbif_sp_scientificName))) cat(cc$lightCoral("Failed to remove NA values. \n")) else cat(cc$lightGreen("NAs removed successfully. \n"))
    } else {
      cat(cc$lightGreen("No NAs found. \n"))
    }
    
    if (any(gbif_sp_scientificName == "")) {
      cat(cc$lightCoral("Some species are blank, removing... \n"))
      gbif_sp_scientificName <- data.table(gbif_sp_scientificName[rowSums(gbif_sp_scientificName != "") > 0, ])
      if (any(is.na(gbif_sp_scientificName))) cat(cc$lightCoral("Failed to remove blank values. \n")) else cat(cc$lightGreen("blank values removed successfully. \n"))
    } else {
      cat(cc$lightGreen("No blank values found. \n"))
    }
    
    cat("Removing duplicates. \n")
    
    # Remove species that are already present in the Arctic
    gbif_sp_scientificName <-  anti_join(gbif_sp_scientificName, arctic_present, by = c("species" = "refinedScientificName"))
    
    # Remove species already a part of the glonaf_absent species
    gbif_sp_scientificName <-  anti_join(gbif_sp_scientificName, arctic_absent, by = c("species" = "refinedScientificName"))
    
    # Run synonym Check on these names, then run the entire process again
    if (file.exists("./outputs/filter/gbif-acquisition/wfo-one-uniq.csv")) {
      cat("GBIF synonym check already conducted, reading file... \n")
      gbif_species <- fread("./outputs/filter/gbif-acquisition/wfo-one-uniq.csv")
    } else {
      names(gbif_sp_scientificName$scientificName) <- "rawName"
      
      gbif_sp_wfo <- check_syn_wfo(
        checklist = gbif_sp_scientificName, 
        column = "rawName", 
        folder = "./outputs/filter/gbif-acquisition", 
        max.cores = cores.max, 
        verbose = T, 
        counter = 100
      )
      
      gbif_sp_wfo_one <- check_syn_wfo_one(
        wfo_checklist = gbif_sp_wfo, 
        folder = "./outputs/filter/gbif-acquisition"
      )
      
      gbif_species <- gbif_sp_wfo_one$wfo_one_uniq
      
    }
    
    gbif_sp_list <- select_wfo_column(
      filepath = "./outputs/filter/gbif-acquisition/wfo-one-uniq.csv", 
      col.unique = c("genus", "specificEpithet", "infraspecificEpithet"), 
      col.select = "scientificNameAuthorship",
      verbose = T
    )
    
    names(gbif_sp_list) <- "gbif_species"
    
    gbif_species <- gbif_sp_list$gbif_species
    
    
    ##################################################
    #           gbif arctic present absent           #
    ##################################################
    
    
    cat("Separating", cc$lightSteelBlue(nrow(gbif_species)), "gbif species into arctic present and absent. \n")
    
    # First merge to only get species from both dfs
    gbif_present <- merge(gbif_species, arctic_present, by = "refinedScientificName")
    
    cat("GBIF has",  cc$lightSteelBlue(nrow(gbif_present)), "species already present in the Arctic. \n")
    
    # Remove species that are already present in the Arctic
    gbif_absent <-  anti_join(gbif_species, arctic_present, by = "refinedScientificName")
    
    # Remove species already a part of the arctic_absent species
    gbif_arc_absent <-  anti_join(gbif_absent, arctic_absent, by = "refinedScientificName")
    
    cat("Removed", cc$lightSteelBlue(nrow(gbif_absent) - nrow(gbif_arc_absent)), "GBIF species already existing in the Arctic absent list. \n")
    
    cat("There are", cc$lightSteelBlue(nrow(gbif_arc_absent)), "Unique GBIF species not in the Arctic absent list. \n")
    
    # Update the absent list
    #arctic_absent = union_dfs(gbif_glo_absent, arctic_absent, verbose = T)
    
    
    ##################################################
    #        Occurrence arctic absent download       #
    ##################################################
    
    arctic_absent_keys <- get_sp_keys(
      sp_names = data.table(refinedScientificName = arctic_absent$refinedScientificName),
      out.dir = "./outputs/filter/arctic",
      verbose = T
    )
    
    arctic_absent_occ <- get_occ_data(
      species_w_keys = arctic_absent_keys,
      file.name = "outputs/filter/arctic/arctic-absent-occ", 
      region = NULL, 
      download.key = "0049125-231120084113126",
      download.doi = "https://doi.org/10.15468/dl.yhk883"
    )
    
    chunk_dir = "./outputs/filter/arctic/chunk"
    sp_w_keys_out <- arctic_absent_keys
    sp_occ_path <- "./outputs/filter/arctic/arctic-absent-occ.csv"
    
    if (!is.null(arctic_absent_occ)) sp_occ <- arctic_absent_occ else sp_occ <- NULL
    
    
    ##################################################
    #        Occurrence gbif absent download         #
    ##################################################
    
    # gbif_absent_keys <- get_sp_keys(
    #   sp_names = gbif_absent,
    #   out.dir = "./outputs/filter/gbif-occurrence",
    #   verbose = T
    #   )
    # 
    # gbif_absent_occ <- get_occ_data(
    #   species_w_keys = gbif_absent_keys,
    #   file.name = "outputs/filter/gbif-occurrence/gbif-absent-occ", 
    #   region = NULL, 
    #   download.key = NULL, 
    #   download.doi = NULL
    # )
    
  }
  ##################################################
  #       Occurrence test species download         #
  ##################################################
  
  if (!is.null(test)) {
    if (test == "small") {
      cat("test_small data sample: \n")
      print(head(dfs$test_small$refinedScientificName, 3))
      
      out_dir <- "./outputs/filter/test/test-small"
      
      test_small <- data.table(refinedScientificName = dfs$test_small$refinedScientificName)
      
      test_small_keys <- get_sp_keys(
        sp_names = test_small,
        out.dir = out_dir,
        verbose = T
      )
      
      test_small_occ <- get_occ_data(
        species_w_keys = test_small_keys,
        file.name = paste0(out_dir, "/test-small-occ"),
        region = NULL,
        download.key = "0056255-231120084113126",
        download.doi = "https://doi.org/10.15468/dl.9v7s8g"
      )
      
      sp_occ_path <- paste0(out_dir, "/test-small-occ.csv")
      
      sp_occ <- test_small_occ
      chunk_dir <- paste0(out_dir, "/chunk")
      
      sp_w_keys_out <- test_small_keys
      
    } else if (test == "big") {
      cat("test_big data sample: \n")
      print(head(dfs$test_big$refinedScientificName, 3))
      
      out_dir <- "./outputs/filter/test/test-big"
      
      test_big <- data.table(refinedScientificName = dfs$test_big$refinedScientificName)
      
      test_big_keys <- get_sp_keys(
        sp_names = test_big,
        out.dir = out_dir,
        verbose = T
      )
      
      test_big_occ <- get_occ_data(
        species_w_keys = test_big_keys,
        file.name = paste0(out_dir, "/test-big-occ"),
        region = NULL,
        download.key = "0048045-231120084113126",
        download.doi = "https://doi.org/10.15468/dl.t9kk65"
      )
      sp_occ_path <- paste0(out_dir, "/test-big-occ.csv")
      
      sp_occ <- test_big_occ
      chunk_dir <- paste0(out_dir, "/chunk")
      sp_w_keys_out <- test_big_keys
    } else {
      cat(cc$lightCoral("Incorrect test input. \n"))
    }
  }
  
  if (is.null(sp_occ)) {
    chunk_file(
      file_path = sp_occ_path,
      chunk.name = "species",
      chunk.column = c("species", "infraspecificEpithet"), 
      chunk.dir = chunk_dir, 
      chunk.size = 1e6,
      #cores.max = cores.max,
      iterations = NULL,
      verbose = verbose
    )
    
  } else {
    chunk_loaded_df(
      df = sp_occ,
      chunk.name = "species",
      chunk.column = c("species", "infraspecificEpithet"),
      chunk.dir = chunk_dir,
      verbose = verbose
    )
  }
# If chunk_loaded_df has a vector input, then "combined" must be used as chunk.column parameter, else the same as chunk_loaded_df. chunk.name has to be the same for all.
  clean_chunks(
    chunk.name = "species",
    chunk.column = "combined",
    chunk.dir = chunk_dir,
    sp_w_keys = sp_w_keys_out,
    verbose = verbose
  )
  
  cat(cc$lightGreen("Filtering sequence completed successfully. \n"))
  
  return(paste0(chunk_dir, "/species"))
}
