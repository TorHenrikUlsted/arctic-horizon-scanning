source_all("./src/setup/components")

setup_raw_data <- function(column, test = NULL, verbose) {
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
  if(verbose) print(names(checklist))
  
  checked_df <- syncheck_dfs(checklist, column, out.dir = "./outputs/setup/wrangle", verbose = verbose)

  cat(cc$lightGreen("raw data setup completed successfully. \n"))
}

setup_sp <- function(test = T, big_test = F) {
  if (test == T) {
    sp_df <- run_test(big_test)
  } else {
    cat(blue("Initiating full run. \n"))
  }
  
  return(sp_df)
} # end setup
