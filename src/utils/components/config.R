#################
#  Projections  #
#################

laea_crs <- crs("+proj=laea +lon_0=0 +lat_0=90 +datum=WGS84")
laea_crs <- edit_crs(laea_crs, "PROJCRS", "Lambert Azimuthal Equal Area")
laea_crs <- edit_crs(laea_crs, "BASEGEOGCRS", "WGS 84")

longlat_crs <- crs("+proj=longlat +datum=WGS84 +ellps=WGS84")
longlat_crs <- edit_crs(longlat_crs, "GEOGCRS", "WGS 84")

mollweide_crs <- crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
mollweide_crs <- edit_crs(mollweide_crs, "PROJCRS", "mollenweide")
mollweide_crs <- edit_crs(mollweide_crs, "BASEGEOGCRS", "WGS 84")

stere_north_crs <- crs("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")
stere_north_crs <- edit_crs(stere_north_crs, "PROJCRS", "stere north")
stere_north_crs <- edit_crs(stere_north_crs, "BASEGEOGCRS", "WGS 84")


#################
#    Climate    #
#################

climate.var = "bio"

climate.res = 2.5

#################
#    Memory     #
#################

mem_total <- get_mem_usage("total")

mem_limit <- mem_total * 0.75

total_cores <- detectCores() * 0.4

#################
#    Species    #
#################

infraEpithet_designations <- c(
  "subsp.", 
  "ssp.", 
  "var.", 
  "f."
)

standard_infraEpithets <- c(
  "ssp." = "subsp.", # the right side is the wanted standard
  "var." = "var.", 
  "f." = "f."
)

ignored_designations <- c(
  "aff.", # For species with affinity to another one (uncertain)
  "agg.",  # aggregated species, no need for this word at the end
  "s. lat.", # sensu lato --> species + lower taxon level
  "coll.", # collection, do not need this word
  "sp." # for species that are uncertain and only have author name after
)

#################
#    ggplot     #
#################

ggplot.filler <- function(gradient = "viridis", scale.variable = "c", limits = NULL, breaks = NULL, labels = NULL, begin = NULL, end = NULL, trans = NULL, guide, na.value = "transparent") {
  tryCatch({
    # Syntax is: "gradient-option"
    split_str <- str_split(gradient, "-")[[1]]
    gradient <- split_str[[1]]
    option <- toupper(split_str[[2]])
    
    args <- list(
      option = option, 
      guide = guide,
      na.value = na.value
    )
    
    if (!is.null(labels)) {
      args$labels <- labels
    }
    
    if (!is.null(limits)) {
      args$limits <- limits
    }
    
    if (!is.null(breaks)) {
      args$breaks <- breaks
    }
    
    if (!is.null(begin)) {
      args$begin <- begin
    }
    
    if (!is.null(end)) {
      args$end <- end
    }
    
    if (!is.null(trans)) {
      args$trans <- trans
    }
    
    if (gradient == "viridis") {
      fun <- paste0("scale_fill_viridis_", scale.variable)
      return(do.call(fun, args))
      
    } else if (gradient == "whitebox") {
      fun <- paste0("scale_fill_whitebox_", scale.variable)
      args$palette <- args$option
      args$option <- NULL
      return(do.call(fun, args))
    }
  }, error = function(e) {
    vebcat("Error when trying to use custom ggplot.filler function.", color = "fatalError")
    stop(e)
  })
}

