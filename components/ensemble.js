/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var cavm = ee.Image("projects/master-thesis-375622/assets/CAVMmap"),
    roi = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-180, -90],
          [180, -90],
          [180, 90],
          [-180, 90],
          [-180, -90]]], null, false);
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var bioVars_df = ee.Image("WORLDCLIM/V1/BIO");


Map.setCenter(-5.386895, 75.386895, 2.0);

//Map.addLayer(cavm, {bands: "b1"});

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

print(cavm)

//Map.addLayer(cavm, "Cavm");
//Map.addLayer(annualMeanTemp, visParams, 'Annual Mean Temperature');
//Map.addLayer(warmestMonth, visParamsWarmestMonth, 'Warmest Month');

