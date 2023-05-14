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

## -------------------- ABA -----------------------------------------------

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
