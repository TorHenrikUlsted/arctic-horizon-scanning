wrangle_glonaf <- function(name, column, verbose = FALSE) {
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)

  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")

  if (file.exists(formatted_out) && file.exists(absent_out)) {
    absent <- fread(absent_out)
  } else {
    vebcat("Initiating", name, "wrangling protocol", color = "funInit")
    preformat <- fread(paste0("./resources/data-raw/", name, ".csv"), header = TRUE)
    
    ## format the dataframe
    formatted <- preformat[, .(taxon_orig, standardized_name, author, family_tpl, region_id, status)]
    
    ## find all different statuses
    vebcat("Statuses found:", unique(formatted$status), veb = verbose)
    
    formatted <- formatted[, .(taxon_orig, standardized_name, author)]
    formatted <- unique(formatted, by = "standardized_name")
    data.table::setnames(formatted, c("taxon_orig", "standardized_name", "author"), c("verbatimName", column, paste0(column,"Authorship")))
    
    fwrite(formatted, formatted_out, row.names = F, bom = T)
    fl <- nrow(formatted)
    absent <- formatted
    
    catn("Unique standardized names:", fl)
    
    catn("Writing out to:", colcat(dir, color = "output"))
    
    fwrite(absent, absent_out, row.names = F, bom = T)
    
    al <- nrow(absent)
    n_dups <- fl - al
    
    md_dt <- data.table(
      formatted = fl,
      absent = al,
      duplicate = n_dups,
      lost = fl - (al + n_dups)
    )
    
    mdwrite(
      config$files$post_seq_md,
      text = "3;Wrangled GloNAF",
      data = md_dt
    )
    
    vebcat(name, "wrangling protocol successfully completed.", color = "funSuccess")
  }
  
  absent[, sourceDataset := name]

  return(list(
    absent = absent
  ))
}
