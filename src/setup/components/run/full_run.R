run_full <- function() {
  cat(blue("Initiating full run. \n"))
  
  run_syn_check <- FALSE
  syn_path <- NULL
  file_name <- ""
  download_path <- ""
  download_key <- NULL
  doi <- NULL
  need_occ <- FALSE
  

  cat("Loading species list. \n")
  
  sp_df <- fread("./outputs/filtering/filtering_process_outputs/filtered_species.csv", sep = "\t")
  
  if (!file.exists("./outputs/setup/full_run/gbif/full_occ.csv")) {
    need_occ <- TRUE
    download_path <- "./outputs/setup/full_run/gbif"
    file_name <- paste0(download_path, "/", "full_occ")
    download_key <- NULL
    doi <- NULL
  } else {
    download_path <- "./outputs/setup/full_run/gbif"
    file_name <- paste0(download_path, "/", "full_occ")
    download_key <- NULL
    doi <- NULL
  }
  
  if (need_occ == T) {
    cat(cc$lightCoral("Species occurrence data not found. \n"))
    
    if(is.null(download_key)) {
      cat(cc$lightCoral("Download key not found. \n"))
      cat("Initiating GBIF occurrence download. \n")
    }
    
    # get codes
    cat("Getting species keys. \n")
    sp_w_keys <- get_sp_keys(sp_df$scientificName)
    
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
    
    # Might need the chunker here
    sp_df <- fread("./outputs/setup/full_run/gbif/full_occ.csv")
    
    return(sp_df)
  }
  
  cat(light_green("test setup completed successfully. \n"))
}