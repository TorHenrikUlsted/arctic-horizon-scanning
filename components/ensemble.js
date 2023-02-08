/**** Start of imports. If edited, may not auto-convert in the playground. ****/
var cavm = ee.FeatureCollection("projects/master-thesis-375622/assets/aga_circumpolar_geobotanical_2003"),
    glonaf = ee.FeatureCollection("projects/master-thesis-375622/assets/257_9_257_2_GloNAF_Shapefile"),
    bioVars = ee.Image("WORLDCLIM/V1/BIO"),
    tif = ee.Image("projects/master-thesis-375622/assets/raster_cavm_v1");
/***** End of imports. If edited, may not auto-convert in the playground. *****/
var bioClip = bioVars.clip(cavm);
var gloClip = glonaf.filterBounds(cavm);

print(cavm.geometry())
print(cavm.geometry().type());
print(tif);
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


// Principal component analysis
var getPrincipalComponents = function(centered, scale, region) {
 var arrays = centered.toArray();

  // compute the covariance of the bands within the region
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(),
    geometry: region,
    scale: scale,
    maxPixels: 1e9
  });
  
  console.log(region);
  console.log(scale);
  console.log(covar);

  //get arrayCovariance results and cast to an array
  //this represents the band-to-band covariance within the region
  var covarArray = ee.Array(covar.get('array'));

  //perform an eigen analysis and slice apart the values and vectors
  var eigens = covarArray.eigen();

  //this is a P-length vector of Eigenvalues
  var eigenValues = eigens.slice(1, 0, 1);
  //This is a PxP matrix with eigenvectors in rows
  var eigenVectors= eigens.slice(1, 1);

  //convert the array image to 2D arrays for matrix computations
  var arrayImage = arrays.toArray(1);

  //Left multiply the image array by the matrix of eigenvectors
  var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage);

  //turn the square roots of the Eigenvalues into a P-band image
   var sdImage = ee.Image(eigenValues.sqrt())
      .arrayProject([0]).arrayFlatten([getNewBandNames('SD')]);

  //Turn the PCs into a P-band image, normalized by SD
  return principalComponents
    //Throw out an unneeded dimension, [[]] -> []
    .arrayProject[0]
    //Make the one band array image a multi-band image, [] -> image
    .arrayFlatten([getNewBandNames('pc')])
    //Normalize the PCs by their SDs
    .divide(sdImage);
  }
  
  getPrincipalComponents(tif, 1920, cavm.geometry());


//Map.addLayer(annualMeanTemp, visParams, 'Annual Mean Temperature');
//Map.addLayer(warmestMonth, visParamsWarmestMonth, 'Warmest Month');
//Map.addLayer(glonaf);
//Map.addLayer(arcticGlonaf)
//Map.addLayer(cavm);