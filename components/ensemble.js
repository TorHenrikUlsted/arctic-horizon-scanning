/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var cavm = ee.FeatureCollection("projects/master-thesis-375622/assets/aga_circumpolar_geobotanical_2003"),
    bioVars = ee.Image("WORLDCLIM/V1/BIO");
/***** End of imports. If edited, may not auto-convert in the playground. *****/



var bioClip = bioVars.clip(cavm)
Map.centerObject(cavm);

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
