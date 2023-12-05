## Get and Set working directory
tryCatch({
  setwd("arctic-horizon-scanning")
}, error = function(e) {
  cat("Working directory: ", getwd(), "\n")
})  

# Load packages
  pkgs = c(
    "rgbif",
    "dplyr",
    "WorldFlora",
    "geodata",
    "terra",
    "sf",
    "sp",
    "data.table",
    "CoordinateCleaner",
    "reticulate",
    "stringdist",
    "stringr",
    "parallel",
    "progressr",
    "crayon",
    "corrplot",
    "hypervolume"
  )
  
  # Check for outdated packages without considering dependencies
  outdated = intersect(old.packages()[, "Package"], pkgs)
  if (length(outdated) > 0) {
    outdated_packages = paste(outdated, collapse = ", ")
    ## Ask the user if they want to update the outdated packages
    update = readline(prompt = paste("Outdated packages found: \n", outdated_packages, "\nDo you want to update the outdated packages? [y/n] "))
    if (tolower(update) == "y") {
      ## Update outdated packages
      install.packages(outdated, ask = F)
    } else {
      cat("Packages not updated.\n")
    }
  } else {
    cat("All specified packages are up to date.\n")
  }
  
  # Install packages if necessary and load them
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg)
    }
  }
  
  
  # Download and remember WFO data if already downloaded, load file
  if (!file.exists("resources/classification.csv")) {
    message("A progress window pops up here, check your taskbar")
    WFO.download(save.dir = "resources", WFO.remember = TRUE)
    WFO_file = "resources/classification.csv"
  } else {
    WFO_file = "resources/classification.csv"
  }
  
  ##Download wwf Ecoregions of the World
  if(!file.exists("resources/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp")) {
    downld_url = "https://files.worldwildlife.org/wwfcmsprod/files/Publication/file/6kcchn7e3u_official_teow.zip"
    fallbck_url = "https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world"
    paper_url = "https://doi.org/10.1641/0006-3568(2001)051[0933:TEOTWA]2.0.CO;2"
    zip_file = "resources/wwfTerrestrialEcoRegions.zip"
    dest_dir = "resources/wwfTerrestrialEcoRegions"
    official_dir = file.path("resources/wwfTerrestrialEcoRegions", "official")
    
    
    message("WWF Terrestrial EcoRegions does not exist. trying to download...")
    
    tryCatch({
      ## Try to download the file from the primary URL, direct download
      download.file(downld_url, zip_file)
      ## Unzip the downloaded file
      unzip(zip_file, exdir = dest_dir)
      ## move the files to the desired directory
      files = list.files(official_dir, full.names = TRUE)
      file.rename(files, file.path(dest_dir, basename(files)))
      ## Clean up
      unlink(official_dir, recursive = TRUE)
      file.remove(zip_file)
    }, error = function(e) {
      message("Could not find the download file. trying to open website")
      tryCatch({
        ## Open the file in your default web browser
        browseURL(fallbck_url)
      }, error = function(e) {
        message("WebBrowser also failed, trying last resort. Opening the scientific paper.")
        browseURL(paper_url)
      })
    })
  }
  
  # Time tracker
  format_elapsed_time = function(start_time, end_time) {
    ### Calculate elapsed time in hours, minutes and seconds
    elapsed_time_hours = as.numeric(difftime(end_time, start_time, units = "hours"))
    elapsed_time_minutes = as.numeric(difftime(end_time, start_time, units = "mins"))
    elapsed_time_seconds = as.numeric(difftime(end_time, start_time, units = "secs"))
    
    ### Format elapsed time
    formatted_elapsed_time = paste(floor(elapsed_time_hours), "h", floor(elapsed_time_minutes) %% 60, "m", round(elapsed_time_seconds) %% 60, "s")
    
    return(formatted_elapsed_time)
  }
  
  #Time estimator
  etc = function(start_time, total_time, current_index, total_count) {
    # Calculate the average time per item
    avg_time = total_time / current_index
    
    # Estimate the remaining time
    remaining_time = avg_time * (total_count - current_index)
    
    # Display the estimated remaining time
    cat("Estimated time remaining: ", round(remaining_time, 1), " seconds\n")
    flush.console()
  }
  
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
  
  ## remove text that are not identical after a certain word
  ## "!" negates the result in order to keep the result in the resulting dataframe instead of removing it
  filter_rows_after_split_text = function(df, col1, col2, split_text) {
    df %>% rowwise() %>% filter(!grepl(split_text, !!sym(col1)) | 
                                  !grepl(split_text, !!sym(col2)) | 
                                  (length(strsplit(!!sym(col1), split_text)[[1]]) > 1 && 
                                     length(strsplit(!!sym(col2), split_text)[[1]]) > 1 && 
                                     strsplit(!!sym(col1), split_text)[[1]][2] == 
                                     strsplit(!!sym(col2), split_text)[[1]][2]))
    
    return(df)
  }
  
  
  ## remove text that have a certain input in the middle
  filter_rows_around_split_text = function(df, col1, col2, split_text) {
    df %>% filter(!grepl(split_text, !!sym(col1)) & !grepl(split_text, !!sym(col2)))
  }
  
  extract_name = function(x, prefixes = c("var", "subsp.")) {
    
    # Create a regular expression pattern to match the specified prefixes
    prefix_pattern = paste0(prefixes, collapse = "|")
    
    # Use str_extract to extract the name after the specified prefixes
    name = str_extract(x, paste0("(?<=", prefix_pattern, ")\\s*\\S+"))
    return(name)
  } 
  
  custom_colors = function() {
    paleTurquoise = make_style("#AFEEEE")
    aquamarine = make_style("#7FFFD4")
    lightSteelBlue = make_style("#B0C4DE")
    
    
    return(
      list(
        paleTurquoise = paleTurquoise,
        aquamarine = aquamarine,
        lightSteelBlue = lightSteelBlue
      )
    )
  }
  cc <<- custom_colors()
