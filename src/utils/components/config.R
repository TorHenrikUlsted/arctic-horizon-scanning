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

# Scale the visualized plots raster into specified meters cell resolution
raster_scale_m <- 1000


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
  "ssp." = "subsp.", # the right side is the standard
  "subsp." = "subsp.",
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

# Taxonomically sorted order taxons
angiosperms <- c("Amborellales", "Nymphaeales", "Austrobaileyales", "Canellales", "Piperales", "Magnoliales", "Laurales", "Chloranthales", "Acorales", "Alismatales", "Petrosaviales", "Dioscoreales", "Pandanales", "Liliales", "Asparagales", "Arecales", "Commelinales", "Zingiberales", "Poales", "Ceratophyllales", "Ranunculales", "Proteales", "Trochodendrales", "Buxales", "Gunnerales", "Dilleniales", "Saxifragales", "Celastrales", "Oxalidales", "Malpighiales", "Vitales", "Zygophyllales", "Fabales", "Rosales", "Fagales", "Cucurbitales", "Geraniales", "Myrtales", "Crossosomatales", "Picramniales", "Huerteales", "Sapindales", "Malvales", "Brassicales", "Berberidopsidales", "Santalales", "Caryophyllales", "Cornales", "Ericales", "Icacinales", "Metteniusales", "Garryales", "Gentianales", "Boraginales", "Vahliales", "Solanales", "Lamiales", "Aquifoliales", "Asterales", "Escalloniales", "Bruniales", "Paracryphiales", "Dipsacales", "Apiales")

gymnosperms <- c(
  "Cycadales",
  "Ginkgoales",
  "Araucariales",
  "Cupressales",
  "Pinales",
  "Ephedrales",
  "Welwitchiales",
  "Gnetales"
)

pteridophytes <- c(
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

#################
#     Files     #
#################

post_seq_nums <- "./outputs/post-process/sequence-numbers.md" # Do not change the object name
create_file_if(post_seq_nums, keep = TRUE)

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

#################
#    ggplot     #
#################

gradient.config <- list(
  vis.gradient = "viridis-B",
  guide = guide_legend(reverse = FALSE, title.position = "top", label.position = "bottom", nrow = 1),
  na.value = "transparent"
) 

theme.config <- theme(
  text = element_text(family = "Linux Libertine"),
  plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 16, face = "bold.italic"),
  axis.text = element_text(size = 16),
  axis.title.x = element_text(size = 20, hjust = 0.5),
  axis.title.y = element_text(size = 20),
  legend.text = element_text(size = 20),
  legend.title = element_text(size = 22, hjust = 0.5),
  legend.position = "bottom",
)
