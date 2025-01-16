wrangle_glonaf <- function(name, column, verbose = FALSE) {
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)

  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")
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
  glonaf_species <- formatted
  
  catn("Unique standardized names:", fl)
  
  catn("Writing out to:", colcat(dir, color = "output"))

  fwrite(glonaf_species, absent_out, row.names = F, bom = T)
  
  al <- nrow(glonaf_species)
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

  return(list(
    absent = glonaf_species
  ))
}
