wrangle_test <- function(test = "small", column, verbose = F) {
  if (test == "small") {
    test_small <- fread("./resources/data-raw/test-small.csv", sep = "\t")
    
    setnames(test_small, old = "scientificName", new = column)
    
    return(list(
      test_small = test_small
    ))
  } else if (test == "big") {
    test_big <- fread("./resources/data-raw/test-big.csv", sep = "\t")
    
    setnames(test_big, old = "scientificName", new = column)
    
    return(list(
      test_big = test_big
    ))
  }
}