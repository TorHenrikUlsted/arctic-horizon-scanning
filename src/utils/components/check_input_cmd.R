# Check for allowed command line inputs
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
        message("The commands '", paste(global_validCommands, collapse = "', '"), "' are not allowed together or combined with ", paste(local_validCommands, collapse = ", "), ". Please enter valid inputs separated by commas. \n list names: ", paste(setdiff(local_validCommands, global_validCommands), collapse = ", "), "\n multi-action: ", paste(global_validCommands, collapse = ", "), "\n")
      } else {
        message("Invalid input(s): ", paste(invalidInputs, collapse = ", "), "\n", "Please enter one or more of the following commands separated by commas: ", paste(local_validCommands, collapse = ", "), "\n or: ",paste(global_validCommands, collapse = ", "))
      }
    }
  }
  return(inputString)
}