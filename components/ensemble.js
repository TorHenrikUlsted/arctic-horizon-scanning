var bioVars_df = ee.Image("WORLDCLIM/V1/BIO");
var cavm = Image("projects/master-thesis-375622/assets/CAVMmap")

var annualMeanTemperature = bioVars_df.select('bio01');
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

Map.setCenter(71.72, 52.48, 3.0);
Map.addLayer(annualMeanTemperature, visParams, 'Annual Mean Temperature');
Map.addLayer(warmestMonth, visParamsWarmestMonth, 'Warmest Month');

