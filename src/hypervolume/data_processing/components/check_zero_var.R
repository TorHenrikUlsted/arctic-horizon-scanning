check_zero_var <- function(env_data, verbose = T) {
  cat("Checking for zero variance values. \n")
  
  # Find names of species with zero variance in all dimensions
  zero_var_species <- lapply(env_data, function(x) all(apply(x, 2, function(y) var(y, na.rm = TRUE) == 0 || is.na(var(y, na.rm = TRUE)))))
  
  # Print names of species with zero variance in all dimensions
  zero_var_names <- names(zero_var_species)[unlist(zero_var_species)]
  cat(yellow("Warning: Species with zero variance values:", zero_var_names, "\n"))
  
  if(!is.null(zero_var_names)) {
    cat("Attempting to remove. \n")
    for (name in zero_var_names) {
      # Identify the rows where all dimensions have the same value
      same_value_rows <- apply(env_data[[name]], 1, function(y) length(unique(y)) == 1)
      
      # Remove these rows
      env_data[[name]] <- env_data[[name]][!same_value_rows, ]
    }
  } else {
    cat("No zer values found. \n")
  }
  
  # Remove observations with zero variance for each species
  zero_var_species <- lapply(env_data, function(x) all(apply(x, 2, function(y) var(y, na.rm = TRUE) == 0 || is.na(var(y, na.rm = TRUE)))))
  zero_var_names <- names(zero_var_species)[unlist(zero_var_species)]
  
  cat(if(is.null(zero_var_names)) green("All zero variance successfully removed.") else red("Failed to remove:"), cc$lightSteelBlue(zero_var_names), "\n")
  
  return(env_data)
}