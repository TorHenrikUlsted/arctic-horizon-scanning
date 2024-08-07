handle_period <- function(dataset, period, verbose = FALSE) {
  vebcat("Handling period", veb = verbose)

  if (is.null(period)) {
    return(invisible(NULL))
  }

  existing_periods <- c("1981-2010", "2011-2040", "2041-2070", "2071-2100")

  if (dataset == "climatologies" && !(period %in% existing_periods)) {
    vebcat("Period does not exist.", color = "fatalError")
    vebprint(existing_periods, text = "Existing periods:")
    stop("Edit the 'period' parameter")
  }

  if (dataset != "climatologies" && grepl("-", period)) {
    catn("Converting period to a sequence")
    parts <- str_split(period, "-")[[1]]
    first_year <- parts[1]
    last_year <- parts[2]
    period <- seq(first_year, last_year)
  } else {
    vebcat("Period not converted", veb = verbose)
  }

  return(period)
}

handle_climatologies <- function(dataset, variable, variable.sub = NULL, model = NULL, ssp = NULL, period, version, verbose = FALSE) {
  file_names <- c()
  vebcat("Handling climatologies filenames", veb = verbose)

  if (is.null(period)) {
    vebcat("Cannot have period = NULL with this command.", color = "fatalError")
    stop("Change 'period' parameter.")
  }
  
  existing_models <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
  existing_ssps <- c("ssp126", "ssp370", "ssp585")

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
      vebprint(existing_models, text = "Existing models:")
      vebprint(existing_ssps, text = "Existing ssps:")
      catn("For more information, check out: ")
      catn("https://chelsa-climate.org/downloads/")
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

handle_chelsa_filenames <- function(dataset, variable, variable.sub = NULL, model = NULL, ssp = NULL, period, version, verbose = FALSE) {
  # Handle the different monthly and daily names
  month_var <- FALSE
  day_var <- FALSE
  month <- 1:12
  day <- 1:366

  version <- paste0("V.", version)

  monthly_vars <- c("clt", "cmi", "hurs", "pet", "pr", "rsds", "sfcWind", "tas", "tasmax", "tasmin", "vpd")

  if (variable %in% monthly_vars) month_var <- TRUE

  if (dataset == "daily" || dataset == "daily_normals") day_var <- TRUE

  file_names <- c()

  if (dataset == "daily_normals") {
    file_names <- c(file_names, sprintf("CHELSA_stot_pj_%s_%s.tif", days, version))
    
  } else if (!(dataset %in% c("climatologies", "input"))) {
    vebcat("Handling filenames", veb = verbose)
    
    if (!month_var && !day_var) {
      file_names <- c(file_names, sprintf("CHELSA_%s_%s_%s.tif", variable, year, version))
    } else if (month_var) {
      if (!day) {
        file_names <- c(file_names, sprintf("CHELSA_%s_%02d_%s_%s.tif", variable, month, year, version))
      } else {
        if (dataset == "daily_normals") {
          file_names <- c(file_names, sprintf("CHELSA_stot_pj_%s_%s.tif", day, version))
        } else {
          file_names <- c(file_names, sprintf("CHELSA_%s_%02d_%02d_%s_%s.tif", variable, day, month, year, version))
        }
      }
    }
  } else if (dataset == "climatologies") {
    file_names <- handle_climatologies(dataset, variable, variable.sub, model, ssp, period, version, verbose)
  } else if (dataset == "input") {
    catn("Handling input filenames")
    file_names <- variable
  }

  return(file_names)
}

chelsa_query <- function(dataset = "climatologies", variable = "bio", variable.sub = NULL, period = "1981-2010", model = NULL, ssp = NULL, version = "2.1", verbose = FALSE) {
  base_url <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL"
  
  if (dataset == "climatologies") {
    updated_url <- file.path(base_url, dataset, period, variable)
  } else if (dataset == "input") {
    updated_url <- file.path(base_url, dataset)
  } else {
    updated_url <- file.path(base_url, dataset, variable)
  }

  catn(
    "Acquiring", highcat(dataset), "data for", highcat(variable), "variable in period", highcat(period), 
    if (!is.null(model) && !is.null(ssp)) {
      paste0("with model ", highcat(model), " and ssp ", highcat(ssp))
    } else {
      ""
    },
    "at version", highcat(version)
  )
  
  period <- handle_period(dataset, period, verbose)

  file_names <- handle_filenames(dataset, variable, variable.sub, model, ssp, period, version, verbose)

  return(list(
    url = updated_url,
    files = file_names
  ))
}
  
chelsa_download <- function(dataset = "climatologies", variable = "bio", variable.sub = NULL, period = "1981-2010", model = NULL, ssp = NULL, dir.out = "./resources/climate/chelsa", version = "2.1", verbose = FALSE) {
  vebcat("Initiating Chelsa climate data download", color = "funInit")
  
  out_dir <- paste0(dir.out, "/", variable, "/")
  create_dir_if(out_dir)
  
  tryCatch({
    result <- chelsa_query(dataset, variable, variable.sub, period, model, ssp, version, verbose)
  }, error = function(e) {
    vebcat("Failed to query Chelsa data", color = "fatalError")
    stop(e$message)
  })
  
  full_urls <- file.path(result$url, result$files)
  out_files <- paste0(out_dir, result$files)
  files_to_download <- !file.exists(out_files)
  
  if (any(files_to_download)) {
    tryCatch({
      mapply(
        download.file,
        url = full_urls[files_to_download],
        destfile = destinations[files_to_download],
        SIMPLIFY = FALSE
      )
    }, error = function(e) {
      vebcat("Chelsa climate data failed to download", color = "fatalError")
      stop(e$message)
    })
  } else {
    catn("All files already downloaded.")
  }

  vebcat("Chelsa climate data download completed successfully", color = "funSuccess")
}
