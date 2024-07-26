wrangle_boreal <- function(name, column, verbose = FALSE) {
  
  dir <- paste0("./outputs/setup/wrangle/", name)
  log_dir <- paste0(dir, "/logs")
create_dir_if(dir, log_dir)
  
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  
  load_wwf_ecoregions()
  
  wwf_region <- load_region(
    region.file = "./resources/region/wwfTerrestrialEcoRegions/wwf_terr_ecos.shp", 
    verbose = verbose
  ) 
  
  boreal_region <- wwf_region[wwf_region$BIOME == 6, ]
  
  check_coords_orientation(boreal_region)
  
  boreal_wkt <- vect_to_wkt(
    boreal_region, 
    out.file = log_dir,
    min.x = TRUE, 
    max.x = TRUE
  )
  
  check_coords_orientation(boreal_wkt)
  
  gbif_sp_list <- get_sp_list(
    taxon = "Tracheophyta", 
    region = boreal_wkt, 
    region.name = "boreal", 
    file.name = paste0(dir, "/boreal_sp_list"),
    download.key = "0037892-231120084113126", 
    download.doi = "https://doi.org/10.15468/dl.882wum"
  )
  
  sp <- filter_gbif_list(gbif_sp_list = gbif_sp_list, column = column)
  
  catn("Writing out to:", colcat(dir, color = "output"))
  
  sp <- set_df_utf8(sp)
  fwrite(sp, absent_out, row.names = F, bom = T)
    
  return(list(
    absent = sp
  ))
}
