filter_test_known <- function(dts, column, verbose = FALSE) {
  tp <- dts$test_known
  
  setnames(tp, old = colnames(tp), new = column)
  
  return(list(
    present = tp
  ))
}

filter_test_small <- function(known.filtered, dts, column, verbose = FALSE) {
  test_dir <- "./outputs/filter/test-small"
  
  create_dir_if(test_dir)
  
  test_present <- known.filtered$present
  
  test_small <- dts[["test_small"]]
  
  setnames(test_small, old = colnames(test_small), new = column)
  
  vebprint(test_small, text = "test_small dt:")
  
  # Remove known tests from the test lists
  test_small <- write_filter_fun(
    file.out = paste0(test_dir, "/test-small-final.csv"),
    spec.in = test_small,
    fun = function() {
      test_small <- anti_union(test_small, test_present, column)
      return(test_small)
    }
  )
  
  return(list(
    spec = test_small,
    dir = test_dir
  ))
}

filter_test_big <- function(known.filtered, dts, column, verbose = FALSE) {
  test_dir <- "./outputs/filter/test-big"
  
  create_dir_if(test_dir)
  
  test_present <- known.filtered$present
  
  test_big <- dts[["test_big"]]
  
  setnames(test_big, old = colnames(test_big), new = column)
  
  vebprint(test_big, text = "test_big dt:")
  
  # Remove known tests from the test lists
  test_big <- write_filter_fun(
    file.out = paste0(test_dir, "/test-big-final.csv"),
    spec.in = test_big,
    fun = function() {
      # First merge to only get species from both dfs
      test_big <- anti_union(test_big, test_present, column)
      
      return(test_big)
    })
  
  return(list(
    spec = test_big,
    dir = test_dir
  ))
}



filter_test_st <- function(test = "small") {
  if (test == "small") {
    catn("test_small data sample:")
    print(head(dfs$test_small[[column]], 3))
    
    out_dir <- "./outputs/filter/test/test-small"
    
    test_small <- data.table(column = dts$test_small[[column]])
    
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
    
    test_big <- data.table(column = dts$test_big[[column]])
    
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
