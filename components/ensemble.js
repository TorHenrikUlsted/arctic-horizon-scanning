/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var geometry = /* color: #0b4a8b */ee.Geometry.Polygon(
        [[[-172.55721050202632, 63.58549270046944],
          [-168.16267925202632, 65.69245622307304],
          [-167.63533550202632, 69.74685275015699],
          [-49.67301310535925, 83.53247337455899],
          [33.998861894640726, 81.98766548768653],
          [128.92073689464073, 77.71612306446181],
          [-174.12613810535925, 70.53199763685778],
          [-169.21977248352056, 66.06312797917059],
          [-172.77934279602056, 64.06312431581487],
          [175.53678328160638, 60.54826472874104],
          [111.02506453160639, 69.75171118950178],
          [72.00162703160639, 65.33398843408088],
          [42.85281726316306, 65.84251243279421],
          [41.79812976316306, 68.88264723825593],
          [22.63797351316306, 70.99132322209053],
          [-13.045620236836939, 66.41160927410957],
          [-24.11983898683694, 65.48029688078083],
          [-49.43233898683694, 55.13386872593973],
          [-67.9661351368622, 56.516054940641375],
          [-80.2708226368622, 55.533797762713945],
          [-94.5091038868622, 57.94302924662764],
          [-117.8880101368622, 65.03891363700203],
          [-127.48654178341201, 68.34455723602377],
          [-143.658416783412, 66.66390634657031],
          [-161.236541783412, 57.246429488345825],
          [-175.299041783412, 60.77263076358991]]]);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var cavm = ee.FeatureCollection("projects/master-thesis-375622/assets/aga_circumpolar_geobotanical_2003"),
    glonaf = ee.FeatureCollection("projects/master-thesis-375622/assets/257_9_257_2_GloNAF_Shapefile"),
    bioVars = ee.Image("WORLDCLIM/V1/BIO"),
    tif = ee.Image("projects/master-thesis-375622/assets/raster_cavm_v1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var bioClip = bioVars.clip(cavm);
var gloClip = glonaf.filterBounds(cavm);

print(cavm.geometry())
print(cavm.geometry().type());
print(tif.geometry());
Map.centerObject(cavm);
Map.addLayer(tif);
//Export.image.toAsset(tif, "cavmMapPolygon", "Earth Engine files")
//Export.table.toAsset(cavm, "cavmMapMultipolygon", "https://drive.google.com/drive/u/1/folders/1wTxIM5QenDNmproIueldtahpIxl9zdkQ", 645951)

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


//var arcticGlonaf = glonaf.clip(cavm);

Map.setCenter(-5, 75, 2);

var annualMeanTemp = bioClip.select('bio01');
var visParams = {
  min: -230.0,
  max: 300.0,
  palette: ['blue', 'purple', 'cyan', 'green', 'yellow', 'red'],
};

var warmestMonth = bioClip.select('bio05');
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