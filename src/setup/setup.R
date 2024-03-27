source_all("./src/setup/components")

check_cpu_speed <- function(df.path, max.cores, sample.size = NULL, verbose, counter) {
  hostname <- system("hostname", intern = T)
  dir_path <- paste0("./outputs/setup/cpu/", hostname)
  
  if (file.exists(paste0(dir_path, "/estimated_time.txt"))) {
    estimated_time <- readLines(paste0(dir_path, "/estimated_time.txt"))
    
    estimated_time <- as.numeric(estimated_time)
  } else {
    df <- fread(df.path, sep = "\t")
    
    sample_size <- min(nrow(df), sample.size)
    
    # Sample a subset of the data
    subset <- df[sample(nrow(df), size = sample_size), ]
    
    start_time <- Sys.time()
    
    result <- check_syn_wfo(subset, column = colnames(subset), folder = dir_path, min(nrow(df), max.cores), verbose = verbose, counter = counter)
    
    end_time <- Sys.time()
    
    # Calculate the time it took to run your function on the subset
    time_taken <- difftime(end_time, start_time, units = "secs")
    
    # Estimate the time it would take to run the function on the full dataset
    estimated_time <- time_taken / sample_size
    
    writeLines(as.character(estimated_time), paste0(dir_path, "/estimated_time.txt"))
  }
  
  time.const <<- as.numeric(estimated_time)
  
  time.setup <<- as.numeric(30)
  
  catn("Time setup (sec):", highcat(time.setup))
  catn("time constant (sec):", highcat(estimated_time))
}

setup_raw_data <- function(column, test = NULL, max.cores, verbose, counter) {
  vebcat("Setting up raw data.", color = "funInit")
  
  if (!is.null(test) && length(test) > 0) {
    test <- wrangle_test(test = test, column, verbose = verbose)
    
    checklist <- test
    
  } else {
    aba <- wrangle_aba(column, verbose = verbose)
    ambio <- wrangle_ambio(column, verbose = verbose)
    glonaf <- wrangle_glonaf(column, verbose = verbose)
    
    checklist <- c(aba, ambio, glonaf)
  }
  
  vebprint(names(checklist), verbose, "dfs added to checklist:")
  
  checked_dfs <- syncheck_dfs(
    checklist, 
    column,
    out.dir = "./outputs/setup/wrangle", 
    max.cores = max.cores, 
    verbose = verbose, 
    counter = counter
  )
  
  vebcat("Combining no-matches.", veb = verbose)
  
  if (all(sapply(checked_dfs, is.null))) {
    catn("All data frames already exist.")
    
  } else {
    combined_df <- data.frame()
    
    # Loop over each list in checked_dfs
    for(i in 1:length(checked_dfs)){
      if (!is.null(checked_dfs[[i]]$wfo_one_nomatch) && nrow(checked_dfs[[i]]$wfo_one_nomatch) > 0) {
        if (names(checked_dfs)[i] != "") {
          catn("Getting nomatches for:", highcat(names(checked_dfs)[i]))
          
          checked_dfs[[i]]$wfo_one_nomatch$dfOrigin <- names(checked_dfs)[i]
          catn("Adding an origin column with:", highcat(names(checked_dfs)[i]))
        }
        
        # rbind the wfo_one_nomatch data frame from each list
        combined_df <- rbind(combined_df, checked_dfs[[i]]$wfo_one_nomatch)
      }
    }
    
    nomatch_combined <- combined_df[!duplicated(combined_df[[column]]), ]
    
    catn("Writing combined no-matches to:", colcat("./outputs/setup/wrangle", color = "output"))
    
    if (nrow(nomatch_combined) > 0) {
      catn("There were", highcat(nrow(nomatch_combined)), "species without matches. \n")
      fwrite(nomatch_combined, "./outputs/setup/wrangle/combined-wfo-nomatch.csv", bom = T)
    } else {
      catn("There were", highcat(0), "species without any matches. \n")
    }
  }
  
  vebcat("raw data setup completed successfully.", color = "funSuccess")
}

setup_region <- function() {
  region_setup_timer <- start_timer("region-setup-timer")
  
  create_dir_if("./resources/region/promice")
  greenland_file <- "./resources/region/promice/basalmelt.nc"
  
  if (!file.exists(greenland_file)) {
    tryCatch(
      {
        vebcat("Could not find greenland shape.", color = "nonFatalError")
        catn("Downloading...")
        download.file("https://dataverse.geus.dk/api/access/datafile/:persistentId?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B", destfile = greenland_file)
      },
      error = function(e) {
        vebcat("Failed to download greenland shape.", color = "nonFatalError")
        catn(e$message)
        
        vebcat("Opening download page for manual download.", color = "indicator")
        
        Sys.sleep(2000)
        
        utils::browseURL("https://dataverse.geus.dk/file.xhtml?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B&version=2.0")
      }
    )
  }
  
  create_dir_if("./resources/region/glims-db/")
  glims_file <- "./resources/region/glims-db/glims_polygons.shp"
  
  if (!file.exists(glims_file)) {
    vebcat("Could not find glims shape.", color = "nonFatalError")
    catn("Downloading...")
    
    glims_url <- "https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0272_GLIMS_v1/NSIDC-0272_glims_db_north_20230725_v01.0.zip"
    
    download_status <- try(download.file(glims_url, destfile = glims_file), silent = TRUE)
    
    if (inherits(download_status, "try-error")) {
      vebcat("Opening download page for manual download.", color = "indicator")
      
      Sys.sleep(2000)
      
      utils::browseURL("https://nsidc.org/data/nsidc-0272/versions/1")
    }
    
    catn("File downloaded, unzipping...")
    unzip(glims_file, exdir = "./resources/region/glims-db/")
    
    catn("cleaning folder.")
    keep_files <- c("glims_polygons.dbf", "glims_polygons.prj", "glims_polygons.shp", "glims_polygons.shx")
    
    all_files <- list.files("./resources/region/glims_db/")
    
    remove_files <- setdiff(all_files, keep_files)
    
    file.remove(paste0("./resources/region/promice", remove_files))
  }
  
  region <- import_regions(
    shapefiles = c(
      cavm = "./resources/region/cavm2003/cavm.shp",
      greenland = "./resources/region/promice/basalmelt.nc",
      glims = "./resources/region/glims-db/glims_polygons.shp"
    ),
    "./outputs/setup/region/logs"
  )
  
  region$cavm <- reproject_region(region$cavm, projection = "longlat", line_issue = T)
  
  catn("Setting up Greenland shape.")
  crs(region$greenland) <- crs(stere_crs)
  
  region$greenland <- reproject_region(region$greenland, projection = "longlat")
  
  if (identical(crs(region$cavm), crs(region$greenland))) {
    vebcat("Cavm and Greenland crs are identical.", color = "proSuccess")
  } else {
    vebcat("Cavm and Greenland crs are not identical.", color = "nonFatalError")  } 
  
  green_shape <- as.polygons(region$greenland[[9]])
  
  if (any(!is.valid(green_shape))) {
    catn("Some geoms of greenland are invalid.")
    catn("Attempting to fix.")
    
    valid_gl <- makeValid(green_shape)
    if (any(!is.valid(valid_gl))) {
      stop("Failed to fix invalid geoms.")
    } else {
      vebcat("Successfully made all geoms valid.", color = "proSuccess")
    } 
    
  } else {
    catn("All greenland geoms are valid.")
    
    valid_gl <- green_shape
  }
  
  catn("Cropping Greenland Ice to cavm extents.")
  green_cavm <- crop(valid_gl, ext(region$cavm))
  
  catn("Erasing Greenland Ice from cavm.")
  cavm_noice <- erase(region$cavm, green_cavm)
  
  catn("Setting up GLIMS database shape.")
  region$glims <- reproject_region(region$glims, projection = "longlat")
  
  if (identical(crs(cavm_noice), crs(region$glims))) {
    vebcat("Cavm and glims crs are identical.", color = "proSuccess")
  } else {
    vebcat("Cavm and glims crs are not identical.", color = "nonFatalError")
  }
  
  catn("Checking for invalid geoms.")
  if (any(!is.valid(region$glims))) {
    catn("Some geoms of glims are invalid.")
    catn("Attempting to fix.")
    
    valid_glims <- makeValid(region$glims)
    if (any(!is.valid(valid_glims))) {
      stop("Failed to fix invalid geoms.")
    } else {
      vebcat("Successfully made all geoms valid.", color = "proSuccess")
    }
  } else {
    catn("All glims geoms are valid.")
  }
  
  catn("Cropping glims to cavm extents.")
  glims_cavm <- crop(valid_glims, ext(cavm_noice))
  
  catn("Erasing glims from cavm.")
  cavm_noice <- erase(cavm_noice, glims_cavm)
  
  create_dir_if("./outputs/setup/region")
  
  out_shape <- "./outputs/setup/region/cavm-noice.shp"
  
  catn("Writing vector to file:", colcat(out_shape, color = "output"))
  writeVector(cavm_noice, out_shape)
  
  end_timer(region_setup_timer)
  
  vebcat("Region setup completed successfully.", color = "funSuccess")
  
  return(cavm_noice)
}

setup_region_hv <- function(biovars_region, out.dir, name, method) {
  vebcat("Setting up region Hypervolume", color = "funInit")
  
  region_filename <- paste0("./outputs/setup/region/", name, "/hypervolume-", method,".rds")
  
  if (file.exists(region_filename)) {
    catn("Region hypervolume already exists, reading file.")
    region_hv <- readRDS(region_filename)
  } else {
    create_dir_if("./outputs/setup/region/logs")
    create_dir_if(paste0("./outputs/setup/region/", name))
    region_log_out <- paste0("./outputs/setup/region/logs/", name, "-", method, "-output.txt")
    
    create_file_if(region_log_out)
    
    region_log_msg <- paste0("./outputs/setup/region/logs/", name, "-", method, "-message.txt")
    create_file_if(region_log_msg)
    
    try(region_log_out <- file(region_log_out, open = "at"))
    try(region_log_msg <- file(region_log_msg, open = "at"))
    sink(region_log_out, type = "output")
    sink(region_log_msg, type = "message")
    
    region_hv_timer <- start_timer("region_hv_timer")
    
    region_hv <- analyze_region_hv(biovars_region, out.dir, name, method = method, verbose = T)
    
    end_timer(region_hv_timer)
    
    sink(type = "message")
    sink(type = "output")
    close(region_log_out)
    close(region_log_msg)
  }
  
  return(region_hv)
}


setup_hv_sequence <- function(min_disk_space, verbose = T) {
  vebcat("Initiating hypervolume sequence setup.", color = "funInit")
  
  sp_list_setup <- list.files("./outputs/filter/test/test-small/chunk/species", full.names = TRUE)
  hv_setup_dir <- "./outputs/setup/hypervolume"
  create_dir_if(hv_setup_dir)
  setup_check <- paste0(hv_setup_dir, "/logs/setup-check.txt")
  low_file <- paste0(hv_setup_dir, "/logs/peak-mem-low.txt")
  high_file <- paste0(hv_setup_dir, "/logs/peak-mem-high.txt")
  
  
  
  if (file.exists(low_file)) {
    catn("Low peak ram setup already run.")
  } else {
    catn("Running low peak ram setup by using", highcat(sp_list_setup[[1]]), "wait time: 3 min.")
    
    create_file_if(setup_check)
    
    ram_control <- start_mem_tracking(file.out = setup_check, stop_file = paste0(hv_setup_dir, "/stop-file.txt"))
    
    parallell_processing(
      spec.list = sp_list_setup[[1]],
      method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
      accuracy = "accurate",
      hv.projection = "laea",
      proj.incl.t = 0.5,
      iterations = 1,
      min.disk.space = min_disk_space,
      hv.dir = hv_setup_dir,
    )
    
    stop_mem_tracking(ram_control, low_file, paste0(hv_setup_dir, "/stop-file.txt"))
  } 
  
  if (file.exists(high_file)) {
    catn("high peak ram setup already run.")
      
  } else {

    catn("Running high peak ram setup using", highcat(sp_list_setup[[2]]), "wait time: 25 min.")
    create_file_if(setup_check)
    
    ram_control <- start_mem_tracking(file.out = setup_check, stop_file = paste0(hv_setup_dir, "/stop-file.txt"))
    
    #Run a hypervolume sequence of sax. opp.

      parallell_processing(
        spec.list = sp_list_setup[[3]], # list of strings
        method = "box", #box approx 13 min, gaussian 1 hours 10 minutes
        accuracy = "accurate",
        hv.projection = "laea",
        proj.incl.t = 0.5,
        iterations = 1,
        min.disk.space = min_disk_space,
        hv.dir = hv_setup_dir,
      )
    
    stop_mem_tracking(ram_control, high_file, paste0(hv_setup_dir, "/stop-file.txt"))
  }
  
  peak_mem_low <- as.numeric(readLines(low_file))
  peak_mem_high <- as.numeric(readLines(high_file))
  
  rm(sp_list_setup)
  invisible(gc())
  
  vebcat("hypervolume sequence setup completed successfully.", color = "funSuccess")
  
  return(list(
    low = peak_mem_low,
    high = peak_mem_high
  ))
}
