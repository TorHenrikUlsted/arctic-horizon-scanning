/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var geometry = /* color: #98ff00 */ee.Geometry.MultiPolygon(
        [[[[-168.58605792125584, 65.5831020151794],
           [-167.714559028272, 66.53757500342593],
           [-166.29348559359704, 71.41287426099811],
           [-132.19192309359704, 71.63569381167471],
           [-86.13723559359704, 82.85326419360878],
           [-30.590360593597044, 84.23185205866294],
           [1.0502644064029365, 80.05788566958363],
           [46.75338940640293, 82.02124625676218],
           [126.90963940640293, 77.3123114461842],
           [176.19515958312545, 75.26085340724254],
           [-169.20925795352082, 66.4541852224327],
           [-169.42938110428463, 66.07793330744445],
           [-171.27508422928463, 65.37340025441358],
           [-171.53875610428463, 64.46064009831132],
           [-174.8038259067472, 63.79269758260093],
           [175.35242409325284, 61.28134352081326],
           [164.80554909325284, 67.55099647970025],
           [131.75867409325284, 69.35749036822794],
           [103.63367409325284, 71.75367165211122],
           [89.74695534325284, 70.08843718560594],
           [85.52820534325284, 67.88424669409784],
           [69.70789284325284, 65.44960910648444],
           [51.250861593252836, 66.59240333022578],
           [46.50476784325284, 65.66780180434888],
           [27.87195534325284, 70.5618736576058],
           [-13.788200906747159, 65.30312898542961],
           [-24.862419656747157, 66.02741302423567]]],
         [[[-167.98112235686133, 66.29311276640492],
           [-170.33219657561133, 65.39402826553663],
           [-171.47477470061133, 64.10072398653516],
           [-173.45231376311133, 63.24322996650213],
           [-174.26466806418756, 59.24176368447853],
           [-159.76271493918756, 57.959790364284345],
           [-156.07130868918756, 60.47728796082907],
           [-158.88380868918756, 63.25195855411343],
           [-158.18068368918756, 65.63882321458584],
           [-140.60255868918756, 67.89087105187669],
           [-135.76857431418756, 67.15164463124806],
           [-128.64943368918756, 69.08310445917051],
           [-107.29201181418756, 61.16287994545547],
           [-96.64071822609053, 61.1376882385025],
           [-95.14657760109053, 57.74494308314008],
           [-78.18368697609053, 55.37258052086337],
           [-74.84384322609053, 58.57944498614987],
           [-69.65829635109053, 58.11825015737313],
           [-63.418061976090534, 56.31016240487885],
           [-39.42392135109051, 59.39451454441878],
           [-25.009858851090517, 66.68323538833613]]]]);
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