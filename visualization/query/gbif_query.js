

// Gets all occurrences within a given polygon 
//
// Vertices of the polygon is specified as '(x1 y1, x2 y2, x3 y3, ...)'
// Polygons must have anticlockwise ordering of points, or it will give unpredictable results
const getOccurrencesInPolygon = async (vertices) => {
  const response = await fetch(`https://api.gbif.org/v1/occurrence/search?geometry=POLYGON(${vertices})`)

  if(response.status !== 200)
    return console.log('ERROR OCCURED:', response)

  const data = await response.json()

  console.log(data)
}

let polygon = "(-180 68,180 68,180 90,-180 90,-180 68)"
getOccurrencesInPolygon(polygon)

