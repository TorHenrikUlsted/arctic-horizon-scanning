filter_test_known <- function(dts, column, verbose = FALSE) {
  test_dir <- "./outputs/filter/test-known"
  create_dir_if(test_dir)
  
  test_present <- dts$test_known
  
  setnames(test_present, old = colnames(test_present), new = column)
  
  fwrite(test_present, paste0(test_dir, "/test-known-present-final.csv"), bom = TRUE)
  
  return(list(
    present = test_present
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
    file.out = paste0(test_dir, "/test-small-absent-final.csv"),
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
    file.out = paste0(test_dir, "/test-big-absent-final.csv"),
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
