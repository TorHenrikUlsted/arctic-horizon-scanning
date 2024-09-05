wrangle_glonaf <- function(name, column, verbose = FALSE) {
  dir <- paste0("./outputs/setup/wrangle/", name)
  create_dir_if(dir)

  formatted_out <- paste0(dir, "/", name, "-formatted.csv")
  absent_out <- paste0(dir, "/", name, "-absent.csv")
  preformat <- fread(paste0("./resources/data-raw/", name, ".csv"), header = TRUE)

  ## format the dataframe
  formatted <- preformat[, .(standardized_name, family_tpl, region_id, status)]

  ## find all different statuses
  vebcat("Statuses found:", unique(formatted$status), veb = verbose)

  ## get a list of all unique names
  glonaf_species <- data.table(standardized_name = unique(formatted$standardized_name))
  n_dups <- nrow(formatted) - nrow(glonaf_species)
  glonaf_species$standardized_name <- trimws(glonaf_species$standardized_name)

  setnames(glonaf_species, old = "standardized_name", new = column)

  catn("Writing out to:", colcat(dir, color = "output"))

  formatted <- set_df_utf8(formatted)
  fwrite(formatted, formatted_out, row.names = F, bom = T)

  glonaf_species <- set_df_utf8(glonaf_species)
  fwrite(glonaf_species, absent_out, row.names = F, bom = T)
  
  fl <- nrow(formatted)
  al <- nrow(glonaf_species)
  
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
