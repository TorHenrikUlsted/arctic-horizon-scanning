source_all("./src/setup/components")

setup_env <- function(test = T, big_test = F) {
  if (test == T) {
    sp_df <- run_test(big_test)
  } else {
    cat(blue("Initiating full run. \n"))
  }
} # end setup
