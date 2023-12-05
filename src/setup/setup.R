source_all("./src/setup/components")

setup_sp <- function(test = T, big_test = F) {
  if (test == T) {
    sp_df <- run_test(big_test)
  } else {
    cat(blue("Initiating full run. \n"))
  }
  
  return(sp_df)
} # end setup
