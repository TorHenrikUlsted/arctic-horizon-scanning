# Define the function
important_biovars <- function(mat, threshold) {
  # Define the biovariables
  temp_vars <- list(year = c("bio_1", "bio_2", "bio_3", "bio_4", "bio_7"),
                    month = c("bio_5", "bio_6"),
                    quarter = c("bio_8", "bio_9", "bio_10", "bio_11"))
  
  precip_vars <- list(year = c("bio_12", "bio_15"),
                      month = c("bio_13", "bio_14"),
                      quarter = c("bio_16", "bio_17", "bio_18", "bio_19"))
  
  # Function to get the most important biovariable
  get_important <- function(group, mat, threshold) {
    # Get the correlation matrix for the group
    group_mat <- mat[group, group]
    
    # Get the sum of the absolute values of the correlations for each biovariable
    sums <- rowSums(abs(group_mat))
    
    # Order the biovariables by the sums
    ordered <- order(sums, decreasing = TRUE)
    
    # Get the most important biovariable
    important <- group[ordered[1]]
    
    # Remove the most important biovariable from the group
    group <- group[-ordered[1]]
    
    # While there are still biovariables in the group
    while (length(group) > 0) {
      # If the correlation of the next biovariable with the important biovariables is less than the threshold
      if (all(abs(mat[important, group[1]]) < threshold)) {
        # Add it to the important biovariables
        important <- c(important, group[1])
      }
      
      # Remove the biovariable from the group
      group <- group[-1]
    }
    
    return(important)
  }
  
  # Define the function to check each biovariable individually with all others
  check_individual <- function(mat, threshold) {
    # Define the biovariables
    bio_vars <- colnames(mat)
    
    # Initialize an empty list to store the results
    individual <- list()
    
    # For each biovariable
    for (bio in bio_vars) {
      # Get the correlations with all other biovariables
      correlations <- mat[bio, bio_vars]
      
      # Remove the correlation of the biovariable with itself
      correlations <- correlations[-which(bio_vars == bio)]
      
      # Order the biovariables by the absolute correlations
      ordered <- order(abs(correlations), decreasing = TRUE)
      
      # Get the biovariables that have a correlation less than the threshold
      important <- bio_vars[ordered][abs(correlations[ordered]) < threshold]
      
      # Add the important biovariables to the results
      individual[[bio]] <- important
    }
    
    # Convert the list to a data frame
    df_ind <- data.frame(matrix(unlist(individual), nrow=length(individual), byrow=T))
    
    # Set the row names of the data frame to the biovariable names
    rownames(df_ind) <- names(individual)
    
    return(individual)
  }
  
  # Get the most important biovariables for each group
  temp_important <- lapply(temp_vars, get_important, mat = mat, threshold = threshold)
  precip_important <- lapply(precip_vars, get_important, mat = mat, threshold = threshold)
  
  result <- data.frame(tempMain = temp_important$year[1],
                       tempYear = paste(temp_important$year, collapse = ", "),
                       tempMonth = paste(temp_important$month, collapse = ", "),
                       tempQuarter = paste(temp_important$quarter, collapse = ", "),
                       precipMain = precip_important$year[1],
                       precipYear = paste(precip_important$year, collapse = ", "),
                       precipMonth = paste(precip_important$month, collapse = ", "),
                       precipQuarter = paste(precip_important$quarter, collapse = ", "))
  
  # Print the names of the variables in order of importance
  cat("The 4 most important biovariables for temperature are:", result$tempMain, "\n")
  cat("The 4 most important biovariables for precipitation are:", result$precipMain, "\n")
  
  fdir <- "./outputs/data_analysis/correlation/"
  
  if (!dir.exists(fdir)) dir.create(fdir, recursive = T)
  
  # Write the results to CSV files
  fwrite(result, paste0(fdir, "important_biovars_all.csv"), bom = T)
  
  # Check each biovariable individually with all others
  ind_res <- check_individual(mat = mat, threshold = threshold)
  
  # Write the results to CSV files
  #fwrite(ind_res, paste0(fdir, "important_biovars_ind.csv"), bom = T)

  for (bio in rownames(ind_res)) {
    cat("The 4 most important biovariables for", bio, "are:", paste(ind_res[bio, 1:4], collapse = ", "), "\n")
  }
  
  return(list(temp = temp_important, precip = precip_important, individual = ind_res))
}
