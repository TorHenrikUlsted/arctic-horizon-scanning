message("Initiating GloNAF wrangling")
# Here we wrangle with the GloNAF data to remove unecessary columns for the end list
## Read file
glonaf_preformatted = read.csv("resources/glonaf.csv", header = T)

## format the dataframe
glonaf_formatted = select(glonaf_preformatted, c(standardized_name, family_tpl, region_id, status))

## find all different statuses
unique(glonaf_formatted$status)

## get a list of all unique names
glonaf_species = data.frame(standardized_name = unique(glonaf_formatted$standardized_name))
glonaf_species$standardized_name = trimws(glonaf_species$standardized_name)

write.csv(glonaf_preformatted, "outputs/origin_lists/glonaf_preformat.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(glonaf_formatted, "outputs/wrangling_outputs/glonaf/glonaf_formatted.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(glonaf_species, "outputs/wrangling_outputs/glonaf/glonaf_species_names.csv", row.names = F, fileEncoding = "UTF-8")

cat("GloNAF data wrangling completed \n")