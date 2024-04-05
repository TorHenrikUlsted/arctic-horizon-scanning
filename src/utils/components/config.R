#################
#  Projections  #
#################

laea_crs <- crs("+proj=laea +lon_0=0 +lat_0=90 +datum=WGS84")

longlat_crs <- crs("+proj=longlat +datum=WGS84 +ellps=WGS84")

mollweide_crs <- crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

stere_crs <- crs("+proj=stere +lon_0=-45 +lat_0=90 +k=1 +R=6378273 +no_defs")

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
