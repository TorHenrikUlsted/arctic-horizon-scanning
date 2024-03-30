filter_test <- function(test = "small") {
  if (test == "small") {
    catn("test_small data sample:")
    print(head(dfs$test_small[[column]], 3))
    
    out_dir <- "./outputs/filter/test/test-small"
    
    test_small <- data.table(column = dfs$test_small[[column]])
    
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
    catn("test_big data sample:")
    print(head(dfs$test_big[[column]], 3))
    
    out_dir <- "./outputs/filter/test/test-big"
    
    test_big <- data.table(column = dfs$test_big[[column]])
    
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
    vebcat("Incorrect test input.", color = "nonFatalError")
  }
  
  return()
}