message("Initiating Test wrangling")

test_preformatted = read.csv("resources/test_species.csv", header = T)
test_formatted = test_preformatted
test_formatted$Species_SubSpecies = trimws(test_preformatted$Species_SubSpecies)
test_species = test_formatted

if (!file.exists("outputs/origin_lists/")) dir.create("outputs/origin_lists/", recursive = T, showWarnings = T)
if (!file.exists("outputs/wrangling_outputs/test")) dir.create("outputs/wrangling_outputs/test", recursive = T, showWarnings = T)

write.csv(test_preformatted, "outputs/origin_lists/test_preformatted.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(test_formatted, "outputs/wrangling_outputs/test/test_formatted.csv", row.names = F, fileEncoding = "UTF-8")
write.csv(test_species, "outputs/wrangling_outputs/test/test_species.csv", row.names = F, fileEncoding = "UTF-8")

cat("Test data wrangling completed \n")