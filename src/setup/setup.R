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
  
  cat("Time setup (sec):", cc$lightSteelBlue(time.setup), "\n")
  cat("time constant (sec):", cc$lightSteelBlue(estimated_time), "\n")
}

setup_raw_data <- function(column, test = NULL, max.cores, verbose, counter) {
  cat(blue("Setting up raw data. \n"))
  
  if (!is.null(test) && length(test) > 0) {
    test <- wrangle_test(test = test, column, verbose = verbose)
    
    checklist <- test
    
  } else {
    aba <- wrangle_aba(column, verbose = verbose)
    ambio <- wrangle_ambio(column, verbose = verbose)
    glonaf <- wrangle_glonaf(column, verbose = verbose)
    
    checklist <- c(aba, ambio, glonaf)
  }
  
  if (verbose) cat("dfs added to checklist: \n")
  if (verbose) print(names(checklist))
  
  checked_dfs <- syncheck_dfs(
    checklist, 
    column,
    out.dir = "./outputs/setup/wrangle", 
    max.cores = max.cores, 
    verbose = verbose, 
    counter = counter
  )
  
  if (verbose) cat("Combining no-matches. \n")
  
  if (all(sapply(checked_dfs, is.null))) {
    cat("All data frames already exist. \n")
    
  } else {
    combined_df <- data.frame()
    
    # Loop over each list in checked_dfs
    for(i in 1:length(checked_dfs)){
      if (!is.null(checked_dfs[[i]]$wfo_one_nomatch) && nrow(checked_dfs[[i]]$wfo_one_nomatch) > 0) {
        if (names(checked_dfs)[i] != "") {
          cat("Getting nomatches for:", cc$lightSteelBlue(names(checked_dfs)[i]), "\n")
          
          checked_dfs[[i]]$wfo_one_nomatch$dfOrigin <- names(checked_dfs)[i]
          cat("Adding an origin column with:", names(checked_dfs)[i], "\n")
        }
        
        # rbind the wfo_one_nomatch data frame from each list
        combined_df <- rbind(combined_df, checked_dfs[[i]]$wfo_one_nomatch)
      }
    }
    
    nomatch_combined <- combined_df[!duplicated(combined_df[[column]]), ]
    
    cat("Writing combined no-matches to:", yellow("./outputs/setup/wrangle"), "\n")
    
    if (nrow(nomatch_combined) > 0) {
      cat("There were", cc$lightSteelBlue(nrow(nomatch_combined)), "species without matches. \n")
      fwrite(nomatch_combined, "./outputs/setup/wrangle/combined-wfo-nomatch.csv", bom = T)
    } else {
      cat("There were", cc$lightSteelBlue(0), "species without any matches. \n")
    }
  }
  
  cat(cc$lightGreen("raw data setup completed successfully. \n"))
}

setup_sp <- function(test = "small", verbose = F) {
  if (tolower(test) == "small" || tolower(test) == "big") {
    sp_df <- run_test(test)
  } else {
    cat(blue("Initiating full run. \n"))
    sp_df <- filter(verbose = verbose)
  }
  
  return(sp_df)
}

setup_region_hv <- function(biovars_region, name, method) {
  region_filename <- paste0("./outputs/setup/region/", name, "/hypervolume-", method,".rds")
  
  if (file.exists(region_filename)) {
    cat("Region hypervolume already exists. \n")
    region_hv <- readRDS(region_filename)
  } else {
    create_dir_if("./outputs/setup/region/logs")
    create_dir_if(paste0("./outputs/setup/region/", name))
    region_log_out <- paste0("./outputs/setup/region/logs/", name, "-", method, "-output.txt")
    
    if (!file.exists(region_log_out)) {
      file.create(region_log_out)
    } else {
      file.remove(region_log_out)
      file.create(region_log_out)
    }
    
    region_log_msg <- paste0("./outputs/setup/region/logs/", name, "-", method, "-message.txt")
    if (!file.exists(region_log_msg)) {
      file.create(region_log_msg)
    } else {
      file.remove(region_log_msg)
      file.create(region_log_msg)
    }
    
    try(region_log_out <- file(region_log_out, open = "at"))
    try(region_log_msg <- file(region_log_msg, open = "at"))
    sink(region_log_out, type = "output")
    sink(region_log_msg, type = "message")
    
    region_hv_timer <- start_timer("region_hv_timer")
    
    region_hv <- analyze_region_hv(biovars_region, name, method = method, verbose = T)
    
    end_timer(region_hv_timer)
    
    sink(type = "message")
    sink(type = "output")
    close(region_log_out)
    close(region_log_msg)
  }
  
  return(region_hv)
}

setup_region <- function() {
  region_setup_timer <- start_timer("region-setup-timer")
  
  create_dir_if("./resources/region/promice")
  greenland_file <- "./resources/region/promice/basalmelt.nc"
  
  if (!file.exists(greenland_file)) {
    tryCatch(
      {
        cat(red("Could not find greenland shape. \n"))
        cat("Downloading... \n")
        download.file("https://dataverse.geus.dk/api/access/datafile/:persistentId?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B", destfile = greenland_file)
      },
      error = function(e) {
        cat(cc$lightCoral("Failed to download greenland shape. \n"))
        cat(e$message, "\n")
        
        cat(cc$aquamarine("Opening download page for manual download. \n"))
        
        Sys.sleep(2000)
        
        browseURL("https://dataverse.geus.dk/file.xhtml?persistentId=doi:10.22008/FK2/PLNUEO/BNXG0B&version=2.0")
      }
    )
  }
  
  create_dir_if("./resources/region/glims-db/")
  glims_file <- "./resources/region/glims-db/glims_polygons.shp"
  
  if (!file.exists(glims_file)) {
    cat(cc$lightCoral("Could not find glims shape. \n"))
    cat("Downloading... \n")
    
    glims_url <- "https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0272_GLIMS_v1/NSIDC-0272_glims_db_north_20230725_v01.0.zip"
    
    download_status <- try(download.file(glims_url, destfile = glims_file), silent = TRUE)
    
    if (inherits(download_status, "try-error")) {
      cat(cc$aquamarine("Opening download page for manual download. \n"))
      
      Sys.sleep(2000)
      
      utils::browseURL("https://nsidc.org/data/nsidc-0272/versions/1")
    }
    
    cat("File downloaded, unzipping... \n")
    unzip(glims_file, exdir = "./resources/region/glims-db/")
    
    cat("cleaning folder. \n")
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
  
  cat("Setting up Greenland shape. \n")
  crs(region$greenland) <- crs(stere_crs)
  
  region$greenland <- reproject_region(region$greenland, projection = "longlat")
  
  if (identical(crs(region$cavm), crs(region$greenland))) cat(green("Cavm and Greenland crs are identical. \n")) else cat(red("Cavm and Greenland crs are not identical. \n"))
  
  green_shape <- as.polygons(region$greenland[[9]])
  
  if (any(!is.valid(green_shape))) {
    cat(red("Some geoms of greenland are invalid. \n"))
    cat("Attempting to fix \n")
    valid_gl <- makeValid(green_shape)
    if (any(!is.valid(valid_gl))) stop(red("Failed to fix invalid geoms. \n")) else cat(green("Successfully made all geoms valid. \n"))
  } else {
    cat(green("All greenland geoms are valid. \n"))
    valid_gl <- green_shape
  }
  
  cat("Cropping Greenland Ice to cavm extents. \n")
  green_cavm <- crop(valid_gl, ext(region$cavm))
  
  cat("Erasing Greenland Ice from cavm. \n")
  cavm_noice <- erase(region$cavm, green_cavm)
  
  cat("Setting up GLIMS database shape. \n")
  region$glims <- reproject_region(region$glims, projection = "longlat")
  
  if (identical(crs(cavm_noice), crs(region$glims))) cat(green("Cavm and glims crs are identical. \n")) else cat(red("Cavm and glims crs are not identical. \n"))
  
  cat("Checking for invalid geoms. \n")
  if (any(!is.valid(region$glims))) {
    cat(red("Some geoms of glims are invalid. \n"))
    cat("Attempting to fix \n")
    valid_glims <- makeValid(region$glims)
    if (any(!is.valid(valid_glims))) stop(red("Failed to fix invalid geoms. \n")) else cat(green("Successfully made all geoms valid. \n"))
  } else {
    cat(green("All glims geoms are valid. \n"))
  }
  
  cat("Cropping glims to cavm extents. \n")
  glims_cavm <- crop(valid_glims, ext(cavm_noice))
  
  cat("Erasing glims from cavm. \n")
  cavm_noice <- erase(cavm_noice, glims_cavm)
  
  create_dir_if("./outputs/setup/region")
  
  writeVector(cavm_noice, "./outputs/setup/region/cavm-noice.shp")
  
  end_timer(region_setup_timer)
  
  cat(cc$lightGreen("Region setup completed successfully. \n"))
  
  return(cavm_noice)
}
