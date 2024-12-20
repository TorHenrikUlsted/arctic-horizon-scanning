config <- read_yaml("./config.yaml")

#------------------------#
####      Memory      ####
#------------------------#

config$memory <- list(
  mem_total = get_mem_usage("total"),
  mem_limit = get_mem_usage("total") * 0.75,
  total_cores = detectCores() * 0.4
)

config$species$angiosperms <- load_apg(rank = config$species$clade_rank)
config$species$pteridophytes <- load_ppg(rank = config$species$clade_rank)
config$species$gymnosperms <- load_gpg(rank = config$species$clade_rank)

#------------------------#
####    Projection    ####
#------------------------#

for (projection in names(config$projection$crs)) {
  config$projection$crs[[projection]] <- crs(config$projection$crs[[projection]])
}
rm(projection)

# laea
config$projection$crs$laea <- edit_crs(config$projection$crs$laea, "PROJCRS", "Lambert Azimuthal Equal Area")
config$projection$crs$laea <- edit_crs(config$projection$crs$laea, "BASEGEOGCRS", "WGS 84")

# Aeqd
config$projection$crs$aeqd <- edit_crs(config$projection$crs$aeqd, "PROJCRS", "Azimuthal Equidistant projection")
config$projection$crs$aeqd <- edit_crs(config$projection$crs$aeqd, "BASEGEOGCRS", "WGS 84")

# longlat
config$projection$crs$longlat <- edit_crs(config$projection$crs$longlat, "GEOGCRS", "WGS 84")

# mollweide
config$projection$crs$mollweide <- edit_crs(config$projection$crs$mollweide, "PROJCRS", "mollenweide")
config$projection$crs$mollweide <- edit_crs(config$projection$crs$mollweide, "BASEGEOGCRS", "WGS 84")

# stere north
config$projection$crs$stere_north <- edit_crs(config$projection$crs$stere_north, "PROJCRS", "stere north")
config$projection$crs$stere_north <- edit_crs(config$projection$crs$stere_north, "BASEGEOGCRS", "WGS 84")

#------------------------#
####     Species      ####
#------------------------#

# Automatically add entries without periods
no_period_infraEpithet <- setdiff(
  sub("\\.$", "", unlist(config$species$standard_infraEpithets)),
  names(config$species$standard_infraEpithets)
)

for (des in no_period_infraEpithet) config$species$standard_infraEpithets[[des]] <- paste0(des, ".")
rm(no_period_infraEpithet, des)

#------------------------#
####      Files       ####
#------------------------#

create_file_if(config$files$post_seq_md, keep = TRUE)

#------------------------#
####    Custom font   ####
#------------------------#

download_if(
  out.file = "./resources/font/linux-libertine/LinLibertine_R.ttf",
  download.file.ext = "zip",
  download.direct = "https://dl.dafont.com/dl/?f=linux_libertine",
  download.page = "https://www.dafont.com/linux-libertine.font",
  verbose = FALSE
)

import_font_if(
  font = "Linux Libertine",
  paths = "./resources/font/linux-libertine",
  pattern = "LinLibertine"
)

suppressMessages(loadfonts(device = "all"))

#------------------------#
####  Config handler  ####
#------------------------#

#Handle all inputs and split into parts by "."
config_handler <- function(...) {
  call_args <- match.call(expand.dots = FALSE)$`...`
  args <- list(...)
  names(args) <- sapply(call_args, deparse)
  
  set_nested <- function(lst, keys, value) {
    if (length(keys) == 1) {
      lst[[keys]] <- value
    } else {
      if (is.null(lst[[keys[1]]])) {
        lst[[keys[1]]] <- list()
      }
      lst[[keys[1]]] <- set_nested(lst[[keys[1]]], keys[-1], value)
    }
    return(lst)
  }
  
  for (original_name in names(args)) {
    name <- original_name
    if (grepl(".", name, fixed = TRUE)) {
      name <- gsub(".", "$", name, fixed = TRUE)
    }
    name_parts <- strsplit(name, "$", fixed = TRUE)[[1]]
    
    config$simulation <<- set_nested(config$simulation, name_parts, args[[original_name]])
  }
}
