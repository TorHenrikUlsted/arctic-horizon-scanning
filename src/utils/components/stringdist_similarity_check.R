# Define a function to find similar strings
similarity_check = function(df, colname1, colname2, method, threshold) {
  # Initialize the result_df variable to NULL
  result_df = NULL
  
  # Use a try-catch block to handle potential errors
  tryCatch({
    # Calculate the pairwise distances between all elements in the two specified columns using the specified distance method
    dist_matrix = as.matrix(stringdistmatrix(as.character(df[[colname1]]),as.character(df[[colname2]]), method = method))
    
    # Find the indices of all rows where there is at least one element that is greater than 0 and less than or equal to the threshold
    to_remove = which(apply(dist_matrix, 1, function(x) any(x > 0 & x <= threshold)))
    
    # Extract the corresponding names from the first specified column
    result = df[[colname1]][to_remove]
    
    # Get the list of similar names for each name in the result vector
    similar_names_list = lapply(to_remove, function(x) {
      similar_indices = which(dist_matrix[x,] > 0 & dist_matrix[x,] <= threshold)
      similar_names = df[[colname2]][similar_indices]
      return(similar_names)
    })
    
    # Create a data frame where each row contains a pair of similar names
    result_df = data.frame(Name = rep(result, sapply(similar_names_list, length)), Similar_Name = unlist(similar_names_list))
    
    # Remove duplicate pairs of similar names
    result_df = unique(t(apply(result_df, 1, function(x) sort(x))))
    colnames(result_df) = c("name", "similarName")
  }, error = function(e) {
    # Print an error message if an error occurs
    message(paste("An error occurred:", e))
    cat("Expected input is: df, column1, column2, method, threshold \n")
    if (nrow(filtered_species) == 0) {
      message("The filtered_species data frame is empty")
    }
    
    # Check if the filtered_species data frame contains a column named "scientificName"
    if (!("scientificName" %in% colnames(filtered_species))) {
      message("The filtered_species data frame does not contain a column named 'scientificName'")
    }
  })
  
  return(result_df)
}