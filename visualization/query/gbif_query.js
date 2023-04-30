var cavm = ee.FeatureCollection("projects/master-thesis-375622/assets/aga_circumpolar_geobotanical_2003")
let search = async (group, searchTerm) => {
  let response = await fetch(`https://api.gbif.org/v1/${group}/search?${searchTerm}`)
  let data = await response.json()

  if (data.status !== 200) {
    console.log(this.response)
  }
  console.log(data)

  
}

let group = 'occurrence'
let searchTerm = 'Pan troglodytes'

search(group, searchTerm)


