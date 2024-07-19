config <- list(
  
  simulation = list(
    seed = 12034
  ),
  
  projection = list(
    laea = crs("+proj=laea +lon_0=0 +lat_0=90 +datum=WGS84"),
    longlat = crs("+proj=longlat +datum=WGS84 +ellps=WGS84"),
    mollweide = crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"),
    stere_north = crs("+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"),
    raster_scale_m = 1000,
    out = "laea"
  ),
  
  climate = list(
    var = "bio",
    res = 2.5
  ),
  
  memory = list(
    mem_total = get_mem_usage("total"),
    mem_limit = get_mem_usage("total") * 0.75,
    total_cores = detectCores() * 0.4
  ),
  
  species = list(
    file_separator = "_", # Used to replace whitespace in species filenames
    
    infraEpithet_designations = c(
      "subsp.", 
      "ssp.", 
      "var.", 
      "f."
    ),
    
    standard_infraEpithets = c(
      "ssp." = "subsp.", # the right side is the standard
      "subsp." = "subsp.",
      "var." = "var.", 
      "f." = "f."
    ),
    
    standard_infraEpithets_taxonRank = c( # used when chunking into species conservatively
      "SPECIES" = NA,
      "SUBSPECIES" = "subsp.",
      "VARIETY" = "var.",
      "FORM" = "f."
    ),
    
    ignored_designations = c(
      "aff.", # For species with affinity to another one (uncertain)
      "agg.",  # aggregated species, no need for this word at the end
      "s. lat.", # sensu lato --> species + lower taxon level
      "coll.", # collection, do not need this word
      "sp." # for species that are uncertain and only have author name after
    ),
    
    # Taxonomically sorted order taxons
    angiosperms = c("Amborellales", "Nymphaeales", "Austrobaileyales", "Canellales", "Piperales", "Magnoliales", "Laurales", "Chloranthales", "Acorales", "Alismatales", "Petrosaviales", "Dioscoreales", "Pandanales", "Liliales", "Asparagales", "Arecales", "Commelinales", "Zingiberales", "Poales", "Ceratophyllales", "Ranunculales", "Proteales", "Trochodendrales", "Buxales", "Gunnerales", "Dilleniales", "Saxifragales", "Celastrales", "Oxalidales", "Malpighiales", "Vitales", "Zygophyllales", "Fabales", "Rosales", "Fagales", "Cucurbitales", "Geraniales", "Myrtales", "Crossosomatales", "Picramniales", "Huerteales", "Sapindales", "Malvales", "Brassicales", "Berberidopsidales", "Santalales", "Caryophyllales", "Cornales", "Ericales", "Icacinales", "Metteniusales", "Garryales", "Gentianales", "Boraginales", "Vahliales", "Solanales", "Lamiales", "Aquifoliales", "Asterales", "Escalloniales", "Bruniales", "Paracryphiales", "Dipsacales", "Apiales"),
    
    gymnosperms = c(
      "Cycadales",
      "Ginkgoales",
      "Araucariales",
      "Cupressales",
      "Pinales",
      "Ephedrales",
      "Welwitchiales",
      "Gnetales"
    ),
    
    pteridophytes = c(
      "Lycopodiales",
      "Isoetales",
      "Selaginellales",
      "Equisetales",
      "Psilotales",
      "Ophioglossales",
      "Marattiales",
      "Osmundales",
      "Hymenophyllales",
      "Gleicheniales",
      "Schizaeales",
      "Salviniales",
      "Cyatheales",
      "Polypodiales"
    )
  ),
  
  files = list(
    post_seq_md = "./outputs/post-process/markdown/sequence-numbers.md"
  ),
  
  ggplot = list(
    gradient = list(
      vis.gradient = "viridis-B",
      guide = guide_legend(reverse = FALSE, title.position = "top", label.position = "bottom", nrow = 1),
      na.value = "transparent"
    ),
    
    theme = theme(
      text = element_text(family = "Linux Libertine"),
      plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 16, face = "bold.italic"),
      axis.text = element_text(size = 16),
      axis.title.x = element_text(size = 20, hjust = 0.5),
      axis.title.y = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 22, hjust = 0.5),
      legend.position = "bottom",
    )
  )
)

################
#  Projection  #
################

laea_crs <- edit_crs(config$projection$laea, "PROJCRS", "Lambert Azimuthal Equal Area")
laea_crs <- edit_crs(config$projection$laea, "BASEGEOGCRS", "WGS 84")

longlat_crs <- edit_crs(config$projection$longlat, "GEOGCRS", "WGS 84")

mollweide_crs <- edit_crs(config$projection$mollweide, "PROJCRS", "mollenweide")
mollweide_crs <- edit_crs(config$projection$mollweide, "BASEGEOGCRS", "WGS 84")

stere_north_crs <- edit_crs(config$projection$stere_north, "PROJCRS", "stere north")
stere_north_crs <- edit_crs(config$projection$stere_north, "BASEGEOGCRS", "WGS 84")

#################
#     Files     #
#################

create_file_if(config$files$post_seq_md, keep = TRUE)

####################
#    custom font   #
####################

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
