library("geodata")
library("terra")
#library("dynRB")
library("rgbif")
library("factoextra")
library("ggplot2")
library("dplyr")
#library("sf")
library("grateful")

cite_packages(citation.style = "APA7th")


#Get bio variables from worldClim
bioVars = worldclim_global(var="bio", 10, "resources", version="2.1")

#Check bioVars properties 
bioVars
names(bioVars)
ext(bioVars)
crs(bioVars, proj = T, describe = T)
units(bioVars)

#annual temperature [[1]]
bioVars[[1]]
plot(bioVars[[1]])

# ------------------------------------ CAVM shape ------------------------------------------------------------

#Write SpatVector file of cavm shape
# Citation:
# CAVM Team. 2003. Circumpolar Arctic Vegetation Map. (1:7,500,000 scale), Conservation of Arctic Flora and Fauna (CAFF) 
# Map No. 1. U.S. Fish and Wildlife Service, Anchorage, Alaska. ISBN: 0-9767525-0-6, ISBN-13: 978-0-9767525-0-9

cavm = vect("resources/Cavm2003/aga_circumpolar_geobotanical_2003.shp")
crs(cavm, proj = T, describe = T)
plot(cavm)

#Change projection to same as bioVars
cavm = terra::project(cavm, "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +no_defs +units=m")
crs(cavm, proj = T, describe = T)
plot(cavm)

##Separate Arctic into regions level 3 by TDWG
tdwgl3 = vect("resources/tdwgLvl3/level3.shp")
crs(tdwgl3, proj = T, describe = T)
plot(tdwgl3)

#Crop Regions into Arctic
cropAR = crop(tdwgl3, cavm)
plot(cropAR)
cavmRegions = mask(cropAR, cavm)
plot(cavmRegions)
plot(cavm, add =T) #Here it seems some information is lost

#Crop WorldClim data to Arctic CAVM
crop = crop(bioVars, cavm)
#plot(crop)
bioVarsMask = mask(crop, cavm)
plot(bioVarsMask)

#Check out min and max values
terra::minmax(bioVarsMask)

#Check properties and make into dataframe
bioVarsMask
bioVarsMask_df = as.data.frame(bioVarsMask)

#Bio variables explanation:
# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (×100)
# BIO4 = Temperature Seasonality (standard deviation ×10
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter

#----------------------------------------------- PCA ---------------------------------------------

## PCA of bioVars in the Arctic
gPca <- prcomp(bioVarsMask_df, center = TRUE, scale. = TRUE)
##summary describes which principal component accounts for most information
#Use summary to find most important PCs according to proportion of variance
summary(gPca)

scores = as.data.frame(gPca$x)
head(scores[1:4])

##rotation describes which variable accounts for the most information
#Use rotation to find what bio Variable describes the most important PCs best
bioVarsLoadings = gPca$rotation
bioVarsLoadings = as.data.frame(bioVarsLoadings)
names(bioVarsLoadings)

#find biggest absolute values in the different components
bioVarsLoadings$PC1
#find the biggest absolute value in PC1
max(abs(bioVarsLoadings$PC1))
#Find the name of biggest absolute value in PC1
rownames(bioVarsLoadings)[which.max(abs(bioVarsLoadings$PC1))]

#Find max abs names and values for all PCs and put it in a dataframe with all values given
bioVarsMaxAbsPC = do.call(cbind, lapply(seq_along(bioVarsLoadings), function(i) {
  x = bioVarsLoadings[[i]]
  max_i = order(abs(x), decreasing = T)[1:19]
  setNames(data.frame(rownames(bioVarsLoadings)[max_i], abs(x[max_i]), abs( x[max_i]) / sum( abs(x) )*100 ), 
           c(colnames(bioVarsLoadings)[i], paste0(colnames(bioVarsLoadings)[i], "_abs (∝)"), paste0(colnames(bioVarsLoadings)[i], "_pct (%)")))
}))

write.csv(bioVarsMaxAbsPC, "outputs/most Important BioVariables.csv", row.names = F)

plot(gPca)

#loadings only
fviz_pca_var(gPca, 
             axes = c(1, 2 ),
             col.var = "contrib",
             gradient.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
             repel = T,
             )


fviz_pca_var(gPca, 
             axes = c(3, 4),
             col.var = "contrib",
             gradient.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
             repel = T,
)

#a test to visualize it simpler
fviz_pca_var(gPca, labelsize = 4, repel = TRUE,
             select.var = list(cos2 = 0.75), col.var = "contrib", gradient.cols = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"))

#Make it 3D

# Once you have it oriented correctly, run this code:
M <- par3d("userMatrix")
dput(M)

# You'll get a structure as output. Then start your R Markdown document with something like this:
library(rgl)
options(rgl.useNULL = TRUE)
M <- structure(...) # Replace the ... with the structure you got as output

# In each code chunk that produces a plot, write code like this:
plot3d(gPca$rotation[,1:19])
text3d(gPca$rotation[,1:3], texts=rownames(gPca$rotation), col="red", cex=0.8)
coords <- NULL
for (i in 1:nrow(gPca$rotation)) {
  coords <- rbind(coords, rbind(c(0,0,0), gPca$rotation[i,1:3]))
}
lines3d(coords, col="red", lwd=1)
par3d(userMatrix = M)
rglwidget()

## ----------------------------------------------- Gbif --------------------------------------------------------------
#The download from GBIF will be changed at a later date in order to include a whole lot more data
library(rgbif)

#Check extent of cavm
plot(cavm)
lines(ext(cavm))

#convert cavm extent to WKT string
ext(cavm)[1]
ext(cavm)[3]
cavmWKT = sprintf("POLYGON ((%f %f, %f %f, %f %f, %f %f, %f %f))", ext(cavm)[1], ext(cavm)[3], ext(cavm)[1], ext(cavm)[4], 
                  ext(cavm)[2], ext(cavm)[4], ext(cavm)[2], ext(cavm)[3], ext(cavm)[1], ext(cavm)[3])
cavmWKT
#Make it counter-clockwise as per GBIF standard.
cavmWKT = "POLYGON((-169.50929 55.79623,172.06954 55.79623,172.06954 83.62742,-169.50929 83.62742,-169.50929 55.79623))"
#Prepare GBIF keys
taxonName = "Tracheophyta"
taxonKey = name_backbone(taxonName)$usageKey
#Test data download
occ_download_prep(pred("taxonKey", taxonKey), 
                  pred("hasGeospatialIssue", FALSE),
                  pred("hasCoordinate", TRUE),
                  pred("geometry", cavmWKT),
                  pred("occurrenceStatus","PRESENT"),
                  pred_in("basisOfRecord", "Occurrence"),
                  format = "SIMPLE_CSV"
                  )

#Download the data 
vascPlants = occ_download(pred("taxonKey", taxonKey), 
                          pred("hasGeospatialIssue", FALSE),
                          pred("hasCoordinate", TRUE),
                          pred("geometry", cavmWKT),
                          pred("occurrenceStatus","PRESENT"),
                          pred_in("basisOfRecord", "Occurrence"),
                          format = "SIMPLE_CSV"
                          )
# check status
occ_download_wait(vascPlants)

#get the download Data and import to create dataframe
vascPlants_file = occ_download_get(vascPlants, path = "resources", overwrite = T)
vascPlants_df = occ_download_import(vascPlants_file, path = "resources" )

#USE THIS IF ALREADY DOWNLOADED
vascPlants_df = occ_download_import(as.download("resources/0083144-230224095556074.zip"))


## test to see which branch includes NA
gbif_test_na_branch = vascPlants_df %>% 
  group_by(kingdom, phylum, class, order, family, genus) %>% 
  summarise(species = species)

apply(gbif_test_na_branch, 2, function(x) 
  any(is.na(x)) | any(is.infinite(x)) 
  )

## check if any lat or long are NA
with(vascPlants_df, any(is.na(decimalLatitude)))
with(vascPlants_df, any(is.na(decimalLongitude)))
with(vascPlants_df, any(decimalLongitude < -169.50929 | decimalLongitude > 172.06954))
with(vascPlants_df, any(decimalLatitude < 55.79623 | decimalLatitude > 83.62742))

#Make into spatial points
vascPlantsLongLat = cbind(vascPlants_df$decimalLongitude, vascPlants_df$decimalLatitude)
sp_occ = vect(vascPlantsLongLat, vascPlants_df, type="points", crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +no_defs +units=m")
sp_occ

#test plot
plot(cavm)
plot(sp_occ, add=T, col="red")

#Crop GBIF points to Arctic CAVM
#take time -> ETC 43 minutes
startTime = Sys.time()

cropGBIF = crop(sp_occ, cavm)
plot(cropGBIF)
sp_occMask = mask(cropGBIF, cavm)
plot(sp_occMask)

#Check endTime
endTime = Sys.time()
#print the time it took to complete the function
message("Cropping GBIF data to the Arctic complete. Elapsed time: ", endTime - startTime)

#plot cropped data
plot(cavm)
plot(sp_occMask, add=T, col="red")

#Make corpped data into dataframe from 41 900 entries to 10 096 entries
sp_occCavm_df = as.data.frame(sp_occMask)

#uniqe number of classes
n_distinct(sp_occCavm_df$class)
#uniqe number of orders
n_distinct(sp_occCavm_df$order)
#uniqe number of families
n_distinct(sp_occCavm_df$family)
#uniqe number of genuses
n_distinct(sp_occCavm_df$genus)
#uniqe number of species
n_distinct(sp_occCavm_df$species)

#check for empty strings or NA
any(unique(sp_occCavm_df$species == ""))
any(is.na(unique(sp_occCavm_df$species)))

#List of species occurrences in the CAVM
cavmSpList = as.data.frame(unique(sp_occCavm_df$species))
cavmSpList = cavmSpList[!(cavmSpList$`unique(sp_occCavm_df$species)` == ""), ]
cavmSpList = as.data.frame(cavmSpList)
cavmSpList = `colnames<-`(cavmSpList, "species")

#Check for empty strings on new dataframe
any(cavmSpList == "")

#write species list to CSV
write.csv(cavmSpList, "outputs/Species in the CAVM.csv", row.names = F)

# ------------------------------------ CBVM shape ------------------------------------------------------------
cbvm = vect("resources/NABoreal/NABoreal.shp")
plot(cbvm)
crs(cbvm)

# ------------------------------------ Cross reference taxa synonyms ------------------------------------------------------------
library(WorldFlora)

#Download and remember data
WFO.download(save.dir = "resources", WFO.remember = T)

# matching also includes synonyms
wfoArcticMatch = WFO.match(spec.data = sp_occMask, WFO.data = WFO.data,
                        Genus = "genus",
                        Species = "species",
                        Infraspecific.rank = "taxonRank", 
                        Infraspecific = "infraspecificEpithet",
                        Authorship = "verbatimScientificNameAuthorship",
                        )
# Check for gender differences
WFO.acceptable.match(wfoArcticMatch)
# Ignore differences in vowels
WFO.acceptable.match(wfoArcticMatch, no.vowels = T)

# Outputs
accepted.cases = WFO.acceptable.match(wfoArcticMatch, no.vowels = T)
wfoArcticMatch.accepted = wfoArcticMatch[accepted.cases == T, ]
wfoArcticMatch.denied = wfoArcticMatch[accepted.cases == F, ]

# ------------------------------------ Cross reference with ABA species list ------------------------------------------------------------

# read CSV file
ABA_preformat = read.csv("resources/ABA_2013.csv", header = F)
ncol(ABA_preformat)

## format the ABA CSV file
# Remove empty columns
ABA_formatted = select(ABA_preformat, -tail(seq_along(ABA_preformat), 3))

# Assign new column names using two rows
colnames(ABA_formatted) = paste(ABA_formatted[3, ], ABA_formatted[5, ])
# Add new names to columns 29 to 42
colnames(ABA_formatted)[29:39] = c("ArcticOccurence", "ArcticEndemicSpeciesAE", "Borderline", "Introduced", "Naturalized", "nonNativeStableCasual", "StableCasual", "NativeCasual", "PAF", "GenusCount", "FamilyCount")
ABA_formatted = ABA_formatted[-c(3,5), ]
# Remove row 1:6
ABA_formatted = ABA_formatted[-c(1:4), ]

# Assume that ABA_formatted is your existing data frame
# and that the first column contains the mixed Family, Genus, and Species text

# Create new columns for Class, Family, Genus, and Species
ABA_formatted$Class = ""
ABA_formatted$Family = ""
ABA_formatted$Genus = ""
ABA_formatted$Species = ""
ABA_formatted$Subspecies = ""

# Initialize variables to keep track of the current classification levels
current_class = ""
current_family = ""
current_genus = ""

for (i in seq_len(nrow(ABA_formatted))) {
  
  # Get the text from the first column
  line = ABA_formatted[i,1]
  
  # Split the line into its components
  components = strsplit(line, " ")[[1]]
  
  # Check the number of components to determine the classification level
  if (length(components) == 1) {
    # Class level
    current_class = components[1]
    current_family = ""
    current_genus = ""
    
    #skip the first row of every class as it has no info other than the class name
    if (ABA_formatted[i,"Class"] != "") {
      ABA_formatted[i,"Class"] = current_class
    }
    
  } else if (length(components) >= 2) {
    # Family, Genus, or Species level
    
    # Check the number of digits in the first component
    num_digits = nchar(gsub("[^0-9]", "", components[1]))
    
    if (num_digits == 2) {
      # Family level
      current_family = components[2]
      current_genus = ""
      
      # Update the values in the data frame
      ABA_formatted[i,"Class"] = current_class
      ABA_formatted[i,"Family"] = current_family
      
    } else if (num_digits == 4) {
      # Genus level
      current_genus = components[2]
      
      # Update the values in the data frame
      ABA_formatted[i,"Class"] = current_class
      ABA_formatted[i,"Family"] = current_family
      ABA_formatted[i,"Genus"] = current_genus
    } else if (num_digits >= 6) {
      # Species level
      
      # Remove "ssp." prefix and any first capital letter followed by a dot from species name
      species_name = gsub("^(ssp\\.\\s)?([A-Z]\\.\\s)?", "", paste(components[2:length(components)], collapse = " "))
      
      # Check if the species name contains a sub-species designation with or without parenthesis
      if (grepl("\\(ssp\\.\\s.*\\)", species_name)) {
        # If so, assign the sub-species to the subSpecies column
        sub_species = gsub("^.*\\(ssp\\.\\s(.*)\\)$", "\\1", species_name)
        ABA_formatted[i,"Subspecies"] = sub_species
        
        # Remove the sub-species designation from the species name
        species_name = gsub("\\(ssp\\.\\s.*\\)", "", species_name)
      } else if (grepl("ssp\\.\\s", species_name)) {
        # If so, assign the sub-species to the subSpecies column
        sub_species = gsub("^(.*)(ssp\\.\\s)(.*)$", "\\3", species_name)
        ABA_formatted[i,"Subspecies"] = sub_species
        
        # Remove the sub-species designation from the species name
        species_name = gsub("^(.*)(\\s+ssp\\.\\s+.*)$", "\\1", species_name)
      }
      
      # Update the values in the data frame
      ABA_formatted[i,"Class"] = current_class
      ABA_formatted[i,"Family"] = current_family
      ABA_formatted[i,"Genus"] = current_genus
      ABA_formatted[i,"Species"] = species_name
    }
  }
}

#Remove column 1 which is now the old information of class, family, genus, species and subSpecies
ABA_formatted = ABA_formatted[, -1]
# Remove empty rows and rows with all columns NA
ABA_formatted = ABA_formatted[!apply(ABA_formatted, 1, function(x) all(is.na(x) | x == "")),]
#Create new column with Genus and species Combined
ABA_formatted$GenusSpecies = with(ABA_formatted, paste (Genus, Species, sep = " "))
# Move the new columns to be the first four columns of the data frame
ABA_formatted = ABA_formatted[,c("Class", "Family", "Genus", "Species", "Subspecies", "GenusSpecies", setdiff(colnames(ABA_formatted), c("Class", "Family", "Genus", "Species", "Subspecies", "GenusSpecies")))]
#Change to lowercase in class and family for cleaner look
ABA_formatted$Class = toupper(substr(ABA_formatted$Class, 1, 1)) %>% paste0(tolower(substr(ABA_formatted$Class, 2, nchar(ABA_formatted$Class))))
ABA_formatted$Family = toupper(substr(ABA_formatted$Family, 1, 1)) %>% paste0(tolower(substr(ABA_formatted$Family, 2, nchar(ABA_formatted$Family))))
# Replace the square symbol with dash in all columns
ABA_formatted[] = lapply(ABA_formatted, function(x) gsub("", "-", x))

## Cross reference WFO list with the ABA list
names(wfoArcticMatch)
names(ABA_formatted)

#compare species with the ones in the ABA list, ish 90% - 100% similar using stringdist library
if(!require(stringdist)) {
  install.packages("stringdist")
  library(stringdist)
}

# set the Levenshtein distance threshold for similarity
ABA_threshold = 1

# check if each combination exists in ABA_formatted
ABA_check = apply(wfoArcticMatch[, "spec.name.species", drop = FALSE], 1, function(x) {
  any(adist(x[1], ABA_formatted$GenusSpecies, partial = TRUE, ignore.case = TRUE) <= ABA_threshold)
})

# get the list of rows from wfoArcticMatch that are present in ABA_formatted
ABA_present = wfoArcticMatch[ABA_check, ]

# get the list of rows from wfoArcticMatch that are absent in ABA_formatted
ABA_absent = wfoArcticMatch[!ABA_check, ]
message("Comparison of dataframes complete")
