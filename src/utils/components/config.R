laea_crs <- crs("+proj=laea +lon_0=0 +lat_0=90 +datum=WGS84")

longlat_crs <- crs("+proj=longlat +datum=WGS84 +ellps=WGS84")

mollweide_crs <- crs("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

stere_crs <- crs("+proj=stere +lon_0=-45 +lat_0=90 +k=1 +R=6378273 +no_defs")

mem_total <- get_mem_usage("total")

mem_limit <- mem_total * 0.75

total_cores <- detectCores() * 0.4