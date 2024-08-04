handle_period <- function(dataset, period) {
  
  catn("Handling period")
  
  if (is.null(period)) return(invisible(NULL))
  
  if (dataset != "climatologies" && grepl("-", period)) {
    catn("Converting period to a sequence")
    parts <- str_split(period, "-")[[1]]
    first_year <- parts[1]
    last_year <-parts[2]
    period <- seq(first_year, last_year)
  } else {
    catn("Period not converted")
  }
  
  return(period)
}

handle_climatologies <- function(dataset, variable, variable.sub = NULL, model = NULL, ssp = NULL, period, version) {
  
  file_names <- c()
  catn("Handling climatologies filenames")
  
  if (is.null(period)) {
    vebcat("Cannot have period = NULL with this command.", color = "fatalError")
    stop("Change 'period' parameter.")
  }
  
  if (period != "2011-2040" && period != "2041-2070" && period != "2071-2100") {
    if (variable != "5km" && variable != "bio") {
      for (month in 1:12) {
        if (variable == "rsds") {
          file_names <- c(file_names, sprintf("CHELSA_%s_%s_%02d_%s.tif", variable, period, month, version))
        }
        file_names <- c(file_names, sprintf("CHELSA_%s_%02d_%s_%s.tif", variable, month, period, version))
      }
      
    } else if (variable == "bio" && is.null(variable.sub)) {
      for (biovar in 1:19) {
        file_names <- c(file_names, sprintf("CHELSA_%s%s_%s_%s.tif", variable, biovar, period, version))
      }
      
    } else if (variable == "bio" && !is.null(variable.sub)) {
      file_names <- sprintf("CHELSA_%s_%s_%s.tif", variable.sub, period, version)
      
    } else if (variable == "5km") {
      file_names <- sprintf("CHELSA_%s_%s_%s.tif", variable.sub, period, version)
    }
    
  } else {
    if (is.null(model) || is.null(ssp)) {
      vebcat("You need to include a model and ssp for this command.", color = "fatalError")
      stop("Edit the input parameters for 'model' and 'ssp'")
    }
    # Using future models
    if (variable == "bio" && is.null(variable.sub)) {
      for (biovar in 1:19) {
        file_names <- c(file_names, sprintf("CHELSA_%s%s_%s_%s_%s_%s.tif", variable, biovar, period, tolower(model), tolower(ssp), version))
      }
    } else if (variable == "bio" && !is.null(variable.sub)) {
      file_names <- sprintf("CHELSA_%s_%s_%s_%s_%s.tif", variable.sub, period, tolower(model), tolower(ssp), version)
    } else {
      for (month in 1:12) {
        file_names <- sprintf("CHELSA_%s_r1i1p1f1_w5e5_%s_%s_%s_%s_norm.tif", tolower(model), tolower(ssp), variable, month, gsub("-", "_", period))
      }
    }
  }
  
 return(file_names)
}

handle_filenames <- function(dataset, variable, variable.sub = NULL, model = NULL, ssp = NULL, period, version) {
  # Handle the different monthly and daily names
  month = FALSE
  day = FALSE
  
  version <- paste0("V.", version)
  
  monthly_vars <- c("clt", "cmi", "hurs", "pet", "pr", "rsds", "sfcWind", "tas", "tasmax", "tasmin", "vpd")
  
  if (variable %in% monthly_vars) month = TRUE
  
  if (dataset == "daily" || dataset == "daily_normals") day = TRUE
  
  file_names <- c()
  
  if (dataset == "daily_normals") {
    for (day in 1:366) {
      file_names <- c(file_names, sprintf("CHELSA_stot_pj_%s_%s.tif", day, version))
    }
  } else if (!(dataset %in% c("climatologies", "input"))) {
    catn("Handling filenames")
    for (year in period) {
      if (!month && !day) {
        file_names <- c(file_names, sprintf("CHELSA_%s_%s_%s.tif", variable, year, version))
      } else if (month) {
        for (month in 1:12) {
          if (!day) {
            file_names <- c(file_names, sprintf("CHELSA_%s_%02d_%s_%s.tif", variable, month, year, version))
          } else {
            for (day in 1:366) {
              if (dataset == "daily_normals") {
                file_names <- c(file_names, sprintf("CHELSA_stot_pj_%s_%s.tif", day, version)) 
              } else {
                file_names <- c(file_names, sprintf("CHELSA_%s_%02d_%02d_%s_%s.tif", variable, day, month, year, version))
              }
            }
          }
        }
      } 
    }
    
  } else if (dataset == "climatologies") {
    file_names <- handle_climatologies(dataset, variable, variable.sub, model, ssp, period, version)
  } else if (dataset == "input") {
    catn("Handling input filenames")
    file_names <- variable
  } 
  
  return(file_names)
}


chelsa_query <- function(dataset = "annual", variable = "swb", variable.sub = NULL, period = "1981-2010", model = NULL, ssp = NULL, version = "2.1", verbose = FALSE) {
  
  base_url <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL"
  
  fixed_url <- file.path(base_url, dataset, variable)
  
  period <- handle_period(dataset, period)
  
  file_names <- handle_filenames(dataset, variable, variable.sub, model, ssp, period, version)
  
  return(list(
    url = fixed_url,
    files = file_names
  ))
}

chelsa_download <- function(dataset = "annual", variable = "swb", variable.sub = NULL, period = "1981-2010", model = NULL, ssp = NULL, dir.out = "./resources/chelsa", version = "2.1", verbose = FALSE) {
  vebcat("Initiating Chelsa climate data download", color = "funInit")
  
  out_dir <- paste0(dir.out, "/", variable, "/")
  
  if (dir.exists(out_dir)) {
    catn("Chelsea directory already exists.")
    return(invisible())
  }
  
  create_dir_if(out_dir)
  
  result <- chelsa_query(dataset, variable, variable.sub, period, model, ssp, version, verbose)
  
  for (file in result$files) {
    full_url <- file.path(result$url, file)
    download.file(full_url, paste0(out_dir, file))
  }
  
  vebcat("Chelsa climate data downloaded successfully", color = "funSuccess")
}
