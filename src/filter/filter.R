source("./src/utils/utils.R")
source("./src/setup/setup.R")
source("./src/filter/components/select_wfo_column.R")
source("./src/filter/components/union_dfs.R")
source("./src/filter/components/get_gbif_species.R")
source("./src/filter/components/convert_vect_wkt.R")
source("./src/filter/components/calc_coords_orientation.R")
source("./src/hypervolume/data_acquisition/components/import_regions.R")
source("./src/hypervolume/data_processing/components/combine/combine_wkt_anticlockwise.R")

filter = function() {
  cat(blue("Initiating filtering sequence. \n"))
  
  filter_timer = start_timer("filter_timer")
  
  ##################################################
  #              Wranlging setup                   #
  ##################################################
  
  setup_raw_data("preScientificName", test = NULL, verbose = T)
  
  ##################################################
  #                     load dfs                   #
  ##################################################
  
  if (verbose) cat("Loading dfs. \n")
  
  lists <- select_wfo_column(read.dir = "./resources/synonym-checked", column = "scientificName", verbose = T)
  
  ##################################################
  #              combine aba ambio                 #
  ##################################################
  
  aba_present = lists$aba_present
  ambio_present  = lists$ambio_present
  
  arctic_present <- union_dfs(aba_present, ambio_present, verbose = T)
  
  aba_absent = lists$aba_absent
  ambio_absent  = lists$ambio_absent
  
  arctic_absent <- union_dfs(aba_absent, ambio_absent, verbose = T)
  
  
  ##################################################
  #         glonaf arctic present absent           #
  ##################################################
  
  glonaf_species <- lists$glonaf_species
  
  # First merge to only get species from both dfs
  glonaf_present <- merge(glonaf_species, arctic_present, by = "scientificName")
  
  cat(cc$lightSteelBlue(nrow(glonaf_species) - nrow(glonaf_present)), "GloNAF species already exist in the Arcitc. \n")
  
  glonaf_absent <-  anti_join(glonaf_species, arctic_present, by = "scientificName")
  
  # Update the absent list
  arctic_absent = union_dfs(glonaf_absent, arctic_absent, verbose = T)
  
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
  
  gbif_species <- get_gbif_species(boreal_wkt, "boreal", taxon = "Tracheophyta", cores.max = 1, verbose = T)
  
  
  ##################################################
  #           gbif arctic present absent           #
  ##################################################
  
  cat("GBIF has", cc$lightSteelBlue(length(gbif_species)), "species. \n")
  
  # Remove species that are already present in the Arctic
  gbif_absent <-  anti_join(gbif_species, arctic_present, by = "scientificName")
  
  # Remove species already a part of the glonaf_absent species
  gbif_glo_absent <-  anti_join(gbif_absent, glonaf_absent, by = "scientificName")
  
  # Run synonym Check on these names, then run the entire process again
  gbif_sp_wfo <- check_syn_wfo(gbif_glo_absent, "scientificName", folder = "./outputs/filter/acquire-gbif-species/wfo")
  
  gbif_sp_wfo_one <- check_syn_wfo_one(gbif_sp_wfo, "scientificName", folder = "./outputs/filter/acquire-gbif-species/wfo")
  
  gbif_species <- gbif_sp_wfo_one$scientificName
  
  # First merge to only get species from both dfs
  gbif_present <- merge(gbif_species, arctic_present, by = "scientificName")
  
  cat(cc$lightSteelBlue(nrow(gbif_species) - nrow(gbif_present)), "GBIF species already exist in the Arcitc. \n")
  
  # Remove species that are already present in the Arctic
  gbif_absent <-  anti_join(gbif_species, arctic_present, by = "scientificName")
  
  # Remove species already a part of the glonaf_absent species
  gbif_glo_absent <-  anti_join(gbif_absent, glonaf_absent, by = "scientificName")
  
  cat(cc$lightSteelBlue(nrow(gbif_absent) - nrow(gbif_glo_absent)), "GBIF species already exist in the GloNAF database. \n")
  
  # Update the absent list
  arctic_absent = union_dfs(gbif_glo_absent, arctic_absent, verbose = T)
  
  print(nrow(gbif_glo_absent))
  
  
  # source("./filtering/components/gbif_occ_retrieval.R")
  # sp_occ_df <- occ_retrieval(gbif_species_check)
  # 
  # source("./filtering/components/gbif_sp_crop.R")
  # sp_occ_cropped <- crop_species(sp_occ_df)
  
  
  
  end_timer(filter_timer)
}
