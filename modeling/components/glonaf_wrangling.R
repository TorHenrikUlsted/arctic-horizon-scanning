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

cat("GloNAF data wrangling completed \n")