tesdt <- function(data, projection, verbose = F) {
  cat(blue("Initiating hypervloume output data check\n"))
  
  data_item_len <- data.table(item_name = character(), item_length = integer())
  raster_list <- list()
  df_list <- list()
  
  if (verbose) cat("Separating raster and data frame objects.\n")
  
  for (i in seq_along(data)) {
    item <- data[[i]]
    
    cat(cc$aquamarine("Splitting item:", cc$lightSteelBlue(i), "/", cc$lightSteelBlue(length(data))), "\n")
    
    if (verbose) print(class(item))
    
    if ("SpatRaster" %in% class(item)) {
      if (verbose) cat("Checking rasters. \n")
      
      if (verbose) cat("Item layers:", cc$lightSteelBlue(terra::nlyr(item)), "\n")
      
      if (verbose) {
        cat("Data table output of item:\n")
        print(data.table(item_name = names(data)[i], item_len = nlyr(item)))
      }
      
      old_len <- nlyr(item)
      data_item_len <- rbind(data_item_len, data.table(item_name = names(data)[i], item_length = nlyr(item)))
      raster_list <- append(raster_list, list(item))
      names(raster_list)[length(raster_list)] <- names(data)[i]
      new_len <- nlyr(item)
      
      if (old_len == new_len) {
        cat(cc$lightGreen("The lengths of", names(data)[i], "are identical.\n"))
      } else {
        cat(cc$lightCoral("The lengths of", names(data)[i], "are not identical.\n"))
      }
    } else if ("RasterLayer" %in% class(item) || "RasterStack" %in% class(item) || "RasterBrick" %in% class(item) || "stars" %in% class(item)) {
      cat(cc$lightCoral("This function uses the Terra package. Modify the function if you want to use another package. \n"))
    }
    
    if ("data.frame" %in% class(item) || "data.table" %in% class(item)) {
      old_len <- nrow(item)
      
      if (verbose) {
        cat("Checking data frame:", cc$lightSteelBlue(names(data)[i]), "\n")
        cat("Old length of:", cc$lightSteelBlue(old_len), "\n")
        cat("Data table output of item:\n")
        print(data.table(item_name = names(data)[i], item_length = nrow(item)))
      }
      
      data_item_len <- rbind(data_item_len, data.table(item_name = names(data)[i], item_length = nrow(item)))
      df_list <- append(df_list, list(item))
      names(df_list)[length(df_list)] <- names(data)[i]
      
      new_len <- nrow(item)
      if (verbose) cat("New length of:", cc$lightSteelBlue(new_len), "\n")
      
      if (old_len == new_len) {
        cat(cc$lightGreen("The lengths of", names(data)[i], "are identical.\n"))
      } else {
        cat(cc$lightCoral("The lengths of", names(data)[i], "are not identical.\n"))
      }
    }
  }
  
  # Check length of df_list and raster_list
  if (verbose) cat("Length of data:", cc$lightSteelBlue(length(data)), "\n")
  if (verbose) cat("Length of raster_list:", cc$lightSteelBlue(length(raster_list)), "\n")
  if (verbose) cat("Length of dataframes / tables:", cc$lightSteelBlue(length(df_list)), "\n")
  
  if (length(data) > length(raster_list) + length(df_list)) {
    cat(cc$lightCoral("Length of data is bigger than the combined length of the other lists. An item in the list has been lost in the splitting.\n"))
  } else if (length(data) < length(raster_list) + length(df_list)) {
    cat(cc$lightCoral("Length of the splitted lists are bigger than the input data.\n"))
  } else {
    cat(cc$lightGreen("Data splitted sucessfully.\n"))
  }

  r_diff_crs <- list()
  r_diff_ex_crs <- list()
  new_raster_list <- list()

  for (i in seq_along(raster_list)) {
    raster <- raster_list[[i]]
    raster_crs <- as.character(crs(raster))
    layer_ids <- 1:terra::nlyr(raster)

    if (verbose) {
      cat("Checking CRS values for all layer in", cc$lightSteelBlue(names(raster_list)[i]), "\n")
      cat("Raster object:\n")
      print(raster)
      cat("Number of layers caught:", length(layer_ids), "\n")
    } 

    # Check if all layers in the stack have the same crs as the first one
    within_crs_match <- all(raster_crs == raster_crs[1])
    # Check if they match the expected raster projection
    expected_crs_match <- all(raster_crs == as.character(projection))

    dt <- NULL

    # Check for issues within the raster stack
    if (within_crs_match) {
      cat(cc$lightGreen("All layers of", names(raster_list)[i], "Have the same CRS.\n"))
    } else {
      cat(cc$lightCoral("Some of the layers in", names(raster_list)[i], "have different CRS.\n"))

      dt <- data.table(layer_id = layer_ids, within_crs_match = within_crs_match, expected_crs = NA)

      incorrect_layer_ids <- dt[within_crs_match == FALSE]$layer_id

      if (length(incorrect_layer_ids) > 0) {
        cat("There are", length(incorrect_layer_ids), "layer(s) with correct crs.\n")
        incorrect_layers <- terra::subset(raster, incorrect_layer_ids)
      } else {
        cat("There are no layers with the correct CRS.\n")
      }

      correct_layer_ids <- dt[within_crs_match == TRUE]$layer_id

      if (length(correct_layer_ids) > 0) {
        cat("There are", length(correct_layer_ids), "layer(s) with correct crs.\n")
        correct_layers <- terra::subset(raster, correct_layer_ids)
      } else {
        cat("There are no layers with the correct CRS.\n")
      }

      reprojected_layers <- project(incorrect_layers, projection)
      
      combined_rasters <- rbind(correct_layers, reprojected_layers)
      
      new_raster_list[[names(raster_list)[i]]] <- reprojected_layers

      new_raster_list <- list(new_raster_list, c(correct_layers, reprojected_layers))
    }

    # Check for issues for raster stack related to the expected output projection
    if (expected_crs_match) {
      cat(cc$lightGreen("All layers of", names(raster_list)[i], "Have the same CRS as expected CRS.\n"))
    } else {
      cat(cc$lightCoral("Some of the layers in", names(raster_list)[i], "have different CRS than expected CRS.\n"))

      if (is.null(dt)) {
        dt <- data.table(layer_id = layer_ids, expected_crs = expected_crs_match)
      } else {
        dt$expected_crs <- data.table(expected_crs_match)
      }
      
      incorrect_layer_ids <- dt[expected_crs == FALSE]$layer_id
      
      if (length(incorrect_layer_ids) > 0) {
        cat("There are", length(incorrect_layer_ids), "layer(s) with correct crs.\n")
        incorrect_layers <- terra::subset(raster, incorrect_layer_ids)
      } else {
        cat("There are no layers with the correct CRS.\n")
      }
      
      if (verbose) {
        cat("Data table sample of expected_crs_match:\n")
        print(head(dt, 5))
        cat("incorrect_layer_ids:\n")
        str(incorrect_layer_ids)
        print(length(incorrect_layer_ids))
      }
      
      if (length(incorrect_layer_ids) > 0) {
        cat("There are", length(incorrect_layer_ids), "layer(s) with incorrect crs.\n")
        
      } else {
        cat("There are no layers with the incorrect CRS.\n")
      }
      
      correct_layer_ids <- dt[expected_crs == TRUE]$layer_id
      
      if (length(correct_layer_ids) > 0) {
        cat("There are", length(correct_layer_ids), "layer(s) with correct crs.\n")
        correct_layers <- terra::subset(raster, correct_layer_ids)
      } else {
        cat("There are no layers with the correct CRS.\n")
      }
      
      # Reproject each layer that are incorrect
      reprojected_layers <- project(incorrect_layers, projection)
      
      combined_rasters <- rbind(correct_layers, reprojected_layers)

      new_raster_list[[names(raster_list)[i]]] <- reprojected_layers
    }

    cell_size <- terra::cellSize(raster_list[i], unit = "km", lyrs = TRUE, mask = TRUE)

    cat("Cell size(km) for", names(raster_list)[i], ": ", cc$lightSteelBlue(cell_size), "\n")
    cat("Resolution of raster:", cc$lightSteelBlue(terra::res(raster)), "\n")
  }

  cat(cc$lightGreen("Successfully checked hypervloume output data\n"))
}
