get_stats_data <- function(cleaned_data, hv.dir, hv.method, verbose = F, warn, err) {
  vebcat("Initializing stats data.", color = "funInit")
  # Define directories
  vis_stats_dir <- "./outputs/visualize/stats"
  
  if (!file.exists("./outputs/visualize/stats/included-species.csv") || !file.exists("./outputs/visualize/stats/excluded-species.csv") ) {
    stats_src <- paste0(hv.dir, "/stats/box-stats.csv")
    sp_stats <- cleaned_data
    
    # Get clean species csv
    excluded_sp <- sp_stats %>% filter(excluded == TRUE)
    
    fwrite(excluded_sp, paste0(vis_stats_dir, "/uniq-excluded-species.csv"), bom = TRUE)
    
    included_sp <- sp_stats %>% filter(excluded == FALSE)
    
    fwrite(included_sp, paste0(vis_stats_dir, "/uniq-included-species.csv"), bom = TRUE)
    
    # Get the region for the species
    glonaf_wfo_one <- fread("./resources/synonym-checked/glonaf-species-wfo-one.csv")
    
    gwo <- glonaf_wfo_one %>%
      select(rawName.ORIG) %>%
      # Get the refined scientificNames as used in the hypervolume method
      mutate(species = apply(glonaf_wfo_one[, c("genus", "specificEpithet", "infraspecificEpithet"), drop = FALSE], 1, function(x) paste(na.omit(x), collapse = " ")))
    
    gwo$species <- trimws(gwo$species)
    # Swap space with line to match the sp-stats
    gwo$species <- gsub(" ", "-", gwo$species)
    
    colnames(gwo) <- c("origName", "species")
    
    # Remove duplicates
    gwo <- gwo[!duplicated(gwo$species)]
    
    # Match the names with the sp_stats
    gwo_stats <- sp_stats %>%
      left_join(gwo, by = "species")
    
    # Get the region id
    glonaf_orig_df <- fread("./resources/data-raw/glonaf.csv")
    
    godf <- glonaf_orig_df %>%
      select(standardized_name, region_id)
    
    colnames(godf) <- c("origName", "regionId")
    
    glonaf_orig_region <- fread("./resources/region/glonaf-region-desc.csv")
    
    gor <- glonaf_orig_region %>% select(region_id, country_ISO, country)
    
    colnames(gor) <- c("regionId", "countryIso", "country")
    
    # Match the gwo_stats and godf --- Some will have multiple region ids
    matched_glonaf <- merge(gwo_stats, godf, by = "origName")
    
   vebcat("Number of rows when appending region ids:", highcat(nrow(matched_glonaf)), veb = verbose)
    
    # Remove duplicates based on combined region id and standardized name
    matched_dups <- nrow(matched_glonaf[duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ])
    vebcat("Number of duplicates after appending region id:", highcat(matched_dups), veb = verbose)
    
    if (matched_dups > 0) {
      matched_glonaf <- matched_glonaf[!duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ]
      
     vebcat("Number of duplicates after removal:", highcat(nrow(matched_glonaf[duplicated(paste(matched_glonaf$standardized_name, matched_glonaf$region_id)), ])), veb = verbose)
      
      vebcat("Number of rows after finishing appending region ids:", highcat(nrow(matched_glonaf)), veb = verbose)
    }
    
    # Append region info
    matched_stats <- merge(matched_glonaf, gor, by = "regionId")
    
    matched_stats <- matched_stats %>%
      select(species, origName, observations, dimensions, samplesPerPoint, randomPoints, excluded, jaccard, sorensen, fracVolumeSpecies, fracVolumeRegion, overlapRegion, includedOverlap, regionId, countryIso, country)
    
    # Seperate into included and excluded species
    
    included_sp <- matched_stats %>% filter(excluded == FALSE)
    
    vebcat("Number of included species:", highcat(nrow(included_sp)), veb = verbose)
    
    create_dir_if("./outputs/visualize/stats")
    
    fwrite(included_sp, paste0(vis_stats_dir, "/included-species.csv"), bom = TRUE)
    
    excluded_sp <- matched_stats %>% filter(excluded == TRUE)
    
    vebcat("Number of included species:", highcat(nrow(excluded_sp)), veb = verbose)
    
    fwrite(excluded_sp, paste0(vis_stats_dir, "/excluded-species.csv"), bom = TRUE)
    
  } else {
    included_sp <- fread(paste0(vis_stats_dir, "/included-species.csv"))
    
    excluded_sp <- fread(paste0(vis_stats_dir, "/excluded-species.csv"))
  }
  
  vebcat("Stats successfully initialised.", color = "funSuccess")
  
  return(included_sp)
}

get_inclusion_cell <- function(spec.filename, region = NULL, verbose = FALSE) {
  
  init <- function(spec.filename, region, verbose) {
    
    sp_rast <- load_sp_rast(spec.filename)
    
    if (!identical(ext(sp_rast), ext(region))) {
      catn("Cropping region extent to match species")
      region <- crop(region, ext(sp_rast))
    }
    
    catn("Extracting raster from region.")
    sp_region <- terra::extract(sp_rast, region, cells = TRUE, na.rm = FALSE)
    sp_region <- data.table(sp_region)
    
    vebprint(head(sp_region, 3), verbose)
    
    catn("Converting region to data table.")
    region_dt <- as.data.frame(region)
    region_dt <- as.data.table(region_dt)
    region_dt[, ID := .I]
    
    vebprint(head(region_dt, 3), verbose)
    
    catn("Merging raster_dt and region_dt.")
    sp_regions_dt <- merge(sp_region, region_dt[, .(ID, FLOREG, country, floristicProvince)], by = "ID")
    sp_regions_dt[, ID := NULL]
    
    vebprint(head(sp_regions_dt, 3), verbose)
    
    # Merge cell regions with the raster cells
    catn("Extracting raster cells.")
    sp_dt <- terra::extract(sp_rast, ext(sp_rast), cells = TRUE)
    sp_dt <- as.data.frame(sp_dt)
    sp_dt <- as.data.table(sp_dt)
    names(sp_dt) <- c("cell", "richness") # Richness because each cell can only hold 0 or 1 for each species.
    
    catn("Merging raster cells with region by cells.")
    cell_regions_dt <- merge(sp_dt, sp_regions_dt[, .(cell, FLOREG, country, floristicProvince)], by = "cell", all = TRUE)
    
    cell_regions_dt <- unique(cell_regions_dt, by = "cell")
    cell_regions_dt <- cell_regions_dt[order(cell_regions_dt$cell)]
    
    return(cell_regions_dt)
  }
  
  execute <- function(spec.filename, region, verbose) {
    
    catn("Initiating data table.")
    sp_name <- basename(dirname(spec.filename[1]))
    vebcat(sp_name, veb = verbose)
    
    sp_rast <- load_sp_rast(spec.filename[1])
    
    res_dt <- extract_ext_to_dt(sp_rast, value = "valueSum")
    vebprint(res_dt, veb = verbose)
    
    for(i in 2:length(spec.filename)) {
      catn("Execution Iteration:", i)
      species <- spec.filename[i]
      
      sp_name <- basename(dirname(species))
      catn(sp_name)
      
      sp_rast <- load_sp_rast(species)
      
      vebcat("Extracting raster cells.", veb = verbose)
      sp_rast_dt <- extract_ext_to_dt(sp_rast, value = "value")
      
      vebprint(head(sp_rast_dt, 3), veb = verbose)
      
      vebcat("Merging with init results.", veb = verbose)
      # Merge to sum with the richness column
      cell_regions_dt <- merge(res_dt, sp_rast_dt, by = "cell", all = TRUE)
      
      vebprint(head(cell_regions_dt, 3), veb = verbose)
      
      vebcat("Calculating richness", veb = verbose)
      # Sum richness with new value, only if both are not NA, else keep richness value
      cell_regions_dt[, valueSum := ifelse(!is.na(valueSum) & !is.na(value), valueSum + value, ifelse(!is.na(valueSum), valueSum, value)), by = cell]
      
      cell_regions_dt[, value := NULL] # Remove the value column after summing to save ram and space
      res_dt <- cell_regions_dt
      
      vebprint(head(res_dt, 3), veb = verbose)
    }
    
    catn("Finished Executing parallel.")
    
    return(res_dt)
  }
  
  # Handle the parallel results will return a list of results
  process <- function(parallel.res, verbose) { # parallel.res is a list of data tables
    dt_region <- parallel.res[[1]]
    dt_list <- parallel.res[[2]]
    
    summed_dt <- dt_region
    
    vebprint(dt_region, veb = verbose)
    vebprint(dt_list, veb = verbose)
    
    for (i in 1:length(dt_list)) {
      cat("\rProcessing list", i, "/", length(dt_list))
      
      dt_cell <- dt_list[[i]]
      
      merged_dt <- merge(summed_dt, dt_cell, by = "cell")
      
      # Sum the richness
      merged_dt[, richness := ifelse(!is.na(richness) & !is.na(valueSum), richness + valueSum, ifelse(!is.na(richness), richness, valueSum)), by = cell]
      
      merged_dt[, valueSum := NULL]
      
      summed_dt <- merged_dt
    }; catn()
    
    
    return(summed_dt)
  }
  
  return(list(
    init = init,
    execute = execute,
    process = process
  ))
}

convert_template_raster <- function(input.values, hv.dir, hv.method, hv.project.method, projection, projection.method = "near", verbose = F) {
  
  vebcat("Converting raster", color = "funInit")
  sp_dirs <- list.dirs(paste0(hv.dir, "/projections/", hv.method))
  # Remove the directory name
  sp_dirs <- sp_dirs[-1]
  
  sp_filename <- paste0(sp_dirs[[1]], "/", hv.project.method,".tif")
  
  template <- terra::rast(sp_filename)
  # Check length
  catn("Length of input / raster")
  catn(length(input.values), "/", ncell(template))
  
  terra::values(template) <- input.values
  
  template <- check_crs(template, projection = projection, projection.method = projection.method)
  
  vebcat("Raster converted successfully", color = "funSuccess")
  
  return(template)
}

get_world_map <- function(projection, scale = "medium", pole = "north") {
  vebcat("Setting up world Map.", color = "funInit")
  
  world_map <- ne_countries(scale = scale, returnclass = "sf")
  
  world_map <- vect(world_map)
  
  if (pole == "north") {
    equator_extent <- ext(-180, 180, 0, 90)
  } else if (pole == "south") {
    equator_extent <- ext(-180, 180, -90, 0)
  } else {
    stop("pole has to be either 'south' or 'north'.")
  }
  
  world_map <- crop(world_map, equator_extent)
  
  world_map <- project(world_map, projection)
  
  vebcat("World Map setup completed successfully", color = "funSuccess")
  
  return(world_map)
}

get_inc_coverage <- function(spec.filename, region = NULL, verbose = FALSE) {
  
  execute <- function(spec.filename, region, verbose) {
    sp_name <- basename(dirname(spec.filename))
    
    sp_rast <- terra::rast(spec.filename)
    
    sp_coverage <- terra::global(sp_rast, fun = sum, na.rm = TRUE)
    
    sp_dt <- data.table(coverage = sp_coverage, species = sp_name, filename = spec.filename)
    
    setnames(sp_dt, "coverage.sum", "coverage")
    
    return(sp_dt)
  }
  
  return(list(
    execute = execute
  ))
}

stack_projections <- function(filenames, projection, projection.method, binary = FALSE, verbose = FALSE) {
  vebcat("Stacking probability rasters", color = "funInit")
  
  # Initialize raster stack with the first layer
  sp_name <- filenames[['species']][[1]]
  sp_filename <- filenames[['filename']][[1]]
  sp_rast <- terra::rast(sp_filename)
  names(sp_rast) <- sp_name
  
  if(!identical(crs(sp_rast, proj = TRUE), crs(projection, proj = TRUE))) {
    catn("Reprojecting to", as.character(crs(projection, proj = TRUE)))
    stack <- terra::project(sp_rast, projection, method = projection.method)
  } else {
    stack <- sp_rast
  }
  
  for (i in 2:nrow(filenames)) {
    sp_name <- filenames[['species']][[i]]
    sp_filename <- filenames[['filename']][[i]]
    
    cat("\rStacking raster", i, "/", nrow(filenames))
    
    sp_rast <- terra::rast(sp_filename)
    names(sp_rast) <- sp_name
    
    if(!identical(crs(sp_rast, proj = TRUE), crs(projection, proj = TRUE))) {
      if (verbose) catn("Reprojecting to", as.character(crs(projection, proj = TRUE)))
      sp_rast <- terra::project(sp_rast, projection, method = projection.method)
    }
    
    stack <- c(stack, sp_rast)
    
  }; catn()
  
  # Convert values to 0 and 1 if binary
  if (binary) stack <- ceiling(stack)
  
  vebcat("Successfully stacked probability rasters", color = "funSuccess")
  
  return(stack)
}

get_prob_max <- function(spec.filename, region = NULL, verbose = FALSE) {
  
  execute <- function(spec.filename, region, verbose) {
    sp_name <- basename(dirname(spec.filename))
    
    sp_rast <- terra::rast(spec.filename)
    sp_max <- where.max(sp_rast)
    sp_max <- as.data.table(sp_max)
    sp_max <- sp_max[value != 0]
    sp_max <- sp_max[, -(1:2)]
    
    vebprint(sp_max, veb = verbose)
    
    sp_dt <- sp_max[, .(maxValue = value, species = sp_name, filename = spec.filename)]
    
    return(sp_dt)
  }
  
  return(list(
    execute = execute
  ))
}

get_region_richness <- function(spec.filename, region, verbose) {
  
  init <- function(spec.filename, region, verbose) {
    
  }
  
  execute <- function(spec.filename, region, verbose) {
    catn("Setting up raster.")
    sp_name <- basename(dirname(spec.filename))
    
    sp_rast <- terra::rast(spec.filename)
    names(sp_rast) <- sp_name
    
    catn("Extracting raster from region.")
    sp_region <- terra::extract(sp_rast, region, na.rm = FALSE)
    sp_region <- data.table(sp_region)
    
    vebprint(head(sp_region, 3), verbose)
    
    catn("Converting region to data table.")
    region_dt <- as.data.frame(region)
    region_dt <- as.data.table(region_dt)
    region_dt[, ID := .I]
    
    vebprint(head(region_dt, 3), verbose)
    
    catn("Merging raster_dt and region_dt.")
    sp_regions_dt <- merge(sp_region, region_dt[, .(ID, FLOREG, country, floristicProvince)], by = "ID")
    sp_regions_dt[, ID := NULL]
    
    vebprint(head(sp_regions_dt, 3), verbose)
    
    # Calculate abundance for each region
    catn("Calculating abundance.")
    sp_regions_dt[, abundance := sum(.SD[[sp_name]], na.rm = TRUE), by = "FLOREG" ]
    sp_regions_dt <- unique(sp_regions_dt, by = "FLOREG")
    
    vebprint(sp_regions_dt, verbose)
    
    sp_regions_dt <- melt(sp_regions_dt, id.vars = c("FLOREG", "country", "floristicProvince", "abundance"), measure.vars = sp_name, variable.name = "species", value.name = "value")
    
    vebprint(head(sp_regions_dt, 5), verbose)
  }
  
  process <- function(spec.filename, region, verbose) {
    dt[, value := ifelse(abundance >= 1, 1, 0)]
  }
  
  return(list(
    init = init,
    execute = execute,
    process = process
  ))
}

append_taxon <- function(dt, dt.species, verbose = FALSE) {
  vebcat("Appending taxon names to data table.", color = "funInit")
  
  vebprint(dt, veb = verbose)
  vebcat("Number of species:", length(unique(dt[[dt.species]])), veb = verbose)
  
  species <- unique(dt[[dt.species]])
  
  species <- gsub("-", " ", species)
  
  catn("Checking GBIF backbone.\n")
  checklist <- name_backbone_checklist(name_data = species)
  
  checklist <- data.table(checklist[, c("kingdom", "phylum", "order", "family", "genus", "canonicalName", "verbatim_name")])
  
  names(checklist)[7] <- "species"
  
  na_sp <- checklist[is.na(checklist$genus),]
  
  if (nrow(na_sp) > 0) {
    vebcat("Species missing classification.", color = "nonFatalError")
    print(na_sp)
  }
  
  checklist[, species := gsub(" ", "-", species)]
  
  vebprint(checklist, veb = verbose)
  vebprint(dt, veb = verbose)
  
  catn("Merging data tables.")
  merged_dt <- merge(dt, checklist, by = "species")
  
  merged_dt <- data.table(merged_dt)
  
  vebcat("Taxon names successfully appended to data table.", color = "funSuccess")
  
  return(merged_dt)
}

calculate_taxon_richness <- function(dt, taxon, verbose = FALSE) {
  catn("Calculating richness for", highcat(taxon))
  
  taxon_col <- dt[[taxon]]
  
  #dt[!is.na(dt$taxon_col),]
  
  dt[, value := ifelse(abundance >= 1, 1, 0)]
  
  dt[, taxonRichness := sum(abundance, na.rm = TRUE), by = .(FLOREG, get(taxon))]
  
  # Calculate total taxon richness
  dt[, totalRichness := sum(taxonRichness, na.rm = TRUE), by = .(FLOREG)]
  
  # Calculate relative richness
  dt[, relativeRichness := taxonRichness / totalRichness]
  
  return(dt)
}

append_src_country <- function(src.dt, target.dt, level, taxon, verbose) {
  vebcat("Appending source regions to target data table", color = "funInit")
  catn("Getting working groups.")
  tdwg2 <- vect("./resources/region/tdwg/level2/level2.shp")
  tdwg2_dt <- convert_spatial_dt(tdwg2)
  tdwg2_dt <- tdwg2_dt[, .(LEVEL2_NAM, LEVEL2_COD)]
  colnames(tdwg2_dt) <- c("lvl2Name", "lvl2Code")
  
  tdwg3 <- vect("./resources/region/tdwg/level3/level3.shp")
  tdwg3_dt <- convert_spatial_dt(tdwg3)
  tdwg3_dt <- tdwg3_dt[, .(LEVEL3_NAM, LEVEL3_COD, LEVEL2_COD)]
  colnames(tdwg3_dt) <- c("lvl3Name", "lvl3Code", "lvl2Code")
  
  merged_tdwg <- merge(tdwg3_dt, tdwg2_dt, by = "lvl2Code")

  tdwg4 <- vect("./resources/region/tdwg/level4/level4.shp")
  tdwg4_dt <- convert_spatial_dt(tdwg4)
  tdwg4_dt <- tdwg4_dt[, .(Level_4_Na, Level3_cod)]
  colnames(tdwg4_dt) <- c("srcCountry", "lvl3Code")
  
  merged_tdwg <- merge(merged_tdwg, tdwg4_dt, by = "lvl3Code")

  catn("Setting up source country.")
  
  src_dt <- copy(src.dt)
  target_dt <- copy(target.dt)
  
  # Get source country
  src_dt <- src_dt[, .(species, srcCountry = country)]
  
  # combine with tdwg
  src_dt <- merge(src_dt, merged_tdwg, by = "srcCountry")
  
  if (taxon == "species") {
    src_dt <- unique(src_dt, by = "species")
  }

  print(src_dt)
  
  catn("Setting up target region.")
  
  target_dt[, taxonSum := uniqueN(.SD[[taxon]][value == 1], na.rm = TRUE), by = FLOREG]
  
  target_dt <- unique(target_dt, by = c(taxon, "FLOREG"))
  
  target_dt <- target_dt[taxonSum > 0]
  
  catn("Merging data tables.")
  sankey_dt <- merge(target_dt, src_dt, by = taxon, allow.cartesian = TRUE, all = TRUE)
  
  sankey_dt <- unique(sankey_dt, by = c("FLOREG", "srcCountry"))
  
  vebcat("Appended source regions successfully", color = "funSuccess")
  
  return(sankey_dt)
}