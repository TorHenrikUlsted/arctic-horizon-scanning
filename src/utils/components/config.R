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
# order
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

#################
#    ggplot     #
#################

ggtheme.config <- theme(
  plot.title = element_text(color = "black", vjust = -0.5, hjust = 0.5, size = 16, face = "bold.italic"),
  axis.text = element_text(size = 12),
  axis.title.x = element_text(size = 15, hjust = 0.5),
  axis.title.y = element_text(size = 15),
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 15, hjust = 0.5),
  legend.position = "bottom",
)
