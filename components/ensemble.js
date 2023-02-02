var bioVars_df = ee.Image("WORLDCLIM/V1/BIO");
var cavm = ee.Image("projects/master-thesis-375622/assets/CAVMmap");

Map.addLayer(cavm, {bands: "b1"});

var annualMeanTemp = bioVars_df.select('bio01');
var visParams = {
  min: -230.0,
  max: 300.0,
  palette: ['blue', 'purple', 'cyan', 'green', 'yellow', 'red'],
};
var warmestMonth = bioVars_df.select('bio05');
var visParamsWarmestMonth = {
  min: -96,
  max: 490,
  palette: ['blue', 'purple', 'cyan', 'green', 'yellow', 'red'],
};

Map.setCenter(74.997832, 5.386895, 3.0);
//Map.addLayer(annualMeanTemp, visParams, 'Annual Mean Temperature');
//Map.addLayer(warmestMonth, visParamsWarmestMonth, 'Warmest Month');

