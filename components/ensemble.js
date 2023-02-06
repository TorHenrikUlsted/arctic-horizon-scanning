/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var cavm = ee.FeatureCollection("projects/master-thesis-375622/assets/aga_circumpolar_geobotanical_2003"),
    bioVars = ee.Image("WORLDCLIM/V1/BIO"),
    glonaf = ee.FeatureCollection("projects/master-thesis-375622/assets/257_9_257_2_GloNAF_Shapefile"),
    cavmImg = ee.Image("projects/master-thesis-375622/assets/CAVMmap");
/***** End of imports. If edited, may not auto-convert in the playground. *****/

var bioClip = bioVars.clip(cavm);
var gloClip = glonaf.filterBounds(cavm.geometry());
Map.centerObject(cavm);

/*
// Extract geometries from you regions 
// for more than one region (type: featureCollection), do something like:
var regionGeom = cavmImg.map(function(f) {
  return f.geometry();
});

// Now map over your study sites and use intersect to clip them on the region(s)
var studySitesClip = glonaf.map(function(f) {
  return f.intersection(regionGeom, 1); //1 refers to the maxError argument
});
*/

//Map.addLayer(gloClip);

//var arcticGlonaf = glonaf.clip(cavm);

Map.setCenter(-5, 75, 2);

var annualMeanTemp = bioClip.select('bio01');
var visParams = {
  min: -230.0,
  max: 300.0,
  palette: ['blue', 'purple', 'cyan', 'green', 'yellow', 'red'],
};

var warmestMonth = bioVars.select('bio05');
var visParamsWarmestMonth = {
  min: -96,
  max: 490,
  palette: ['blue', 'purple', 'cyan', 'green', 'yellow', 'red'],
};

//Map.addLayer(annualMeanTemp, visParams, 'Annual Mean Temperature');
//Map.addLayer(warmestMonth, visParamsWarmestMonth, 'Warmest Month');
//Map.addLayer(glonaf);
//Map.addLayer(arcticGlonaf)
//Map.addLayer(cavm);