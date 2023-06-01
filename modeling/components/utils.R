# Load packages
## Add packages
packages = c(
  "rgbif",
  "dplyr",
  "WorldFlora",
  "terra",
  "sf",
  "stringdist"
)
## Install packages if necessary and load them
for (package in packages) {
  if (!require(package, character.only = T)) {
    install.packages(package)
    library(package, character.only = T)
  }
}

## Download and remember WFO data if already downloaded, load file
if (!file.exists("resources/classification.csv")) {
  WFO.download(save.dir = "resources", WFO.remember = TRUE)
  WFO_file = "resources/classification.csv"
} else {
  WFO_file = "resources/classification.csv"
}

## Time tracker
format_elapsed_time = function(start_time, end_time) {
  ### Calculate elapsed time in hours and minutes
  elapsed_time_hours = as.numeric(difftime(end_time, start_time, units = "hours"))
  elapsed_time_minutes = as.numeric(difftime(end_time, start_time, units = "mins"))
  
  ### Format elapsed time
  formatted_elapsed_time = paste(round(elapsed_time_hours, 2), "hours (", round(elapsed_time_minutes, 2), "minutes)")
  
  return(formatted_elapsed_time)
}

## Define a function to find similar strings
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
    result_df <- unique(t(apply(result_df, 1, function(x) sort(x))))
    colnames(result_df) <- c("Name", "Similar Name")
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


## Check for allowed command line inputs
checkInputCommands = function(inputString, local_validCommands, global_validCommands = c(), promptMessage) {
  ## Use a while loop to keep the user in the loop until a command is written correctly
  while (length(unlist(lapply(strsplit(inputString, ","), trimws))) == 0 || 
         !all(unlist(lapply(strsplit(inputString, ","), trimws)) %in% c(local_validCommands, global_validCommands)) || 
         (any(global_validCommands %in% unlist(lapply(strsplit(inputString, ","), trimws))) && 
          length(unlist(lapply(strsplit(inputString, ","), trimws))) > 1)) {
    ## commandCheck for which lists to run checks on
    inputString = readline(promptMessage)
    ## make them lowercase
    inputString = tolower(inputString)
    ## Error message
    if (!all(unlist(lapply(strsplit(inputString, ","), trimws)) %in% c(local_validCommands, global_validCommands)) || 
        (any(global_validCommands %in% unlist(lapply(strsplit(inputString, ","), trimws))) && 
         length(unlist(lapply(strsplit(inputString, ","), trimws))) > 1)) {
      ## Print debugging information
      cat("Debugging: \n")
      cat("Inputs are: ", inputString, "\n", "Stored as: ", "\n")
      print(unlist(lapply(strsplit(inputString, ","), trimws)))
      cat("while loop condition:", !any(unlist(lapply(strsplit(inputString, ","), trimws)) %in% local_validCommands), "\n")
      invalidInputs = unlist(lapply(strsplit(inputString, ","), trimws))[!unlist(lapply(strsplit(inputString, ","), trimws)) %in% c(local_validCommands, global_validCommands)]
      if (any(global_validCommands %in% unlist(lapply(strsplit(inputString, ","), trimws))) && length(unlist(lapply(strsplit(inputString, ","), trimws))) > 1) {
        message("The commands '", paste(global_validCommands, collapse = "', '"), "' are not allowed together or combined with ", paste(local_validCommands, collapse = ", "), ". Please enter valid inputs separated by commas. \n local_validCommands: ", paste(setdiff(local_validCommands, global_validCommands), collapse = ", "), "\n global_validCommands: ", paste(global_validCommands, collapse = ", "), "\n")
      } else {
        message("Invalid input(s): ", paste(invalidInputs, collapse = ", "), "\n", "Please enter one or more of the following commands separated by commas: ", paste(local_validCommands, collapse = ", "), "\n or: ",paste(global_validCommands, collapse = ", "))
      }
    }
  }
  return(inputString)
}