wrangle_glonaf <- function(column, verbose = F) {
  cat(blue("Initiating GloNAF wrangling protocol \n"))
  
  glonaf_preformatted = fread("./resources/data-raw/glonaf.csv", header = T)
  
  ## format the dataframe
  glonaf_formatted = select(glonaf_preformatted, c(standardized_name, family_tpl, region_id, status))
  
  ## find all different statuses
  if (verbose) cat("Statuses found:", unique(glonaf_formatted$status), "\n")
  
  ## get a list of all unique names
  glonaf_species = data.frame(standardized_name = unique(glonaf_formatted$standardized_name))
  glonaf_species$standardized_name = trimws(glonaf_species$standardized_name)
  
  setnames(glonaf_species, old = "standardized_name", new = column)
  
  if(verbose) cat("Writing out dfs. \n")
  create_dir_if("./outputs/setup/wrangle/glonaf")
  
  glonaf_formatted <- set_df_utf8(glonaf_formatted)
  fwrite(glonaf_formatted, "./outputs/setup/wrangle/glonaf/glonaf-formatted.csv", row.names = F, bom = T)
  
  create_dir_if("./outputs/filter/wrangle")
  glonaf_species <- set_df_utf8(glonaf_species)
  fwrite(glonaf_species, "./outputs/setup/wrangle/glonaf/glonaf-species.csv", row.names = F, bom = T)
  
  cat(cc$lightGreen("GloNAF wrangling protocol completed successfully. \n"))
 
   return(list(
     glonaf_species = glonaf_species
     ))
}

