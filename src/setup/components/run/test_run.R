run_test <- function(big_test) {
  cat("Initiating test run. \n")

  run_syn_check <- FALSE
  syn_path <- NULL
  file_name <- ""
  download_path <- ""
  download_key <- NULL
  doi <- NULL
  need_occ <- FALSE

  if (big_test == T) {
    cat("Loading big test list. \n")
    
    test_sp <- fread("./resources/test/data_raw/big_test_species.csv", sep = "\t")

    if (!file.exists("outputs/setup/test/gbif/big/wfo_one_checklist.csv")) {
      cat("synonym check needs to be run. \n")

      syn_path <- "./outputs/setup/test/gbif/big"

      run_syn_check <- TRUE
    } else {
      cat("Synonym check already conducted. \n")
      checked_test_sp <- test_sp
    }

    if (!file.exists("./outputs/setup/test/gbif/big/big_test_occ.csv")) {
      need_occ <- TRUE
      download_path <- "./outputs/setup/test/gbif/big"
      file_name <- paste0(download_path, "/", "big_test_occ")
      download_key <- "0003982-231120084113126"
      doi <- "https://doi.org/10.15468/dl.upsbr2"
    } else {
      download_path <- "./outputs/setup/test/gbif/big"
      file_name <- paste0(download_path, "/", "big_test_occ")
      download_key <- "0003982-231120084113126"
      doi <- "https://doi.org/10.15468/dl.upsbr2"
    }
    # ----End of big test ----
  } else if (big_test == F) {
    cat("Loading small test list. \n")
    
    test_sp <- fread("./resources/test/data_raw/small_test_species.csv", sep = "\t")
    
    # Do a synonym_check on them if not already conducted
    if (!file.exists("outputs/setup/test/gbif/small/wfo_one_checklist.csv")) {
      cat("synonym check needs to be run. \n")

      syn_path <- "./outputs/setup/test/gbif/small"

      run_syn_check <- TRUE
    } else {
      cat("Synonym check already conducted. \n")
      checked_test_sp <- test_sp
    }

    if (!file.exists("./outputs/setup/test/gbif/small/small_test_occ.csv")) {
      need_occ <- TRUE
      download_path <- "./outputs/setup/test/gbif/small"
      file_name <- "./outputs/setup/test/gbif/small/small_test_occ"
      download_key <- "0001144-231120084113126"
      doi <- "https://doi.org/10.15468/dl.cmepat"
    } else {
      download_path <- "./outputs/setup/test/gbif/small"
      file_name <- "./outputs/setup/test/gbif/small/small_test_occ"
      download_key <- "0001144-231120084113126"
      doi <- "https://doi.org/10.15468/dl.cmepat"
    }
  } # ----End of small test----

  if (run_syn_check == T) {
    cat("Synonym Check will be saved in:", syn_path, "\n")
    print(names(test_sp))
    names(test_sp)[1] <- "validScientificName"
    
    matched_test_sp <- check_syn_wfo(test_sp, "validScientificName", syn_path)

    checked_test_sp <- check_syn_wfo_one(matched_test_sp, "validScientificName", syn_path)

    run_syn_check <- FALSE
  }

  if (need_occ == T) {
    cat(cc$lightCoral("Species occurrence data not found. \n"))
    
    if(is.null(download_key)) {
      
    }

    # get codes
    cat("Getting species keys. \n")
    sp_w_keys <- get_sp_keys(checked_test_sp$scientificName)

    # Check if the file already exists
    cat("Using keys to download occurrence data. \n")
    sp_df <- get_occ_data(
      sp_w_keys,
      download_path,
      file_name,
      download_key = download_key,
      doi = doi
    )

    need_occ <- FALSE

    return(sp_df)
  } else {
    cat("Species data frame found. \n")
    
    if (big_test == T) sp_df <- fread("./outputs/setup/test/gbif/big/big_test_occ.csv") else  sp_df <- fread("./outputs/setup/test/gbif/small/small_test_occ.csv")

    return(sp_df)
  }

  cat(light_green("test setup completed successfully. \n"))
  
  return(sp_df)
}
