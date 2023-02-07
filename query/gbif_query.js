const search = async (group, searchTerm) => {
  const response = await fetch(`https://api.gbif.org/v1/${group}/search?${searchTerm}`)
  const data = await response.json()

  if (data.status !== 200) {
    console.log(this.response)
  }
  console.log(data)

  
}

var group = 'occurrence'
var searchTerm = 'Pan troglodytes'

search(group, searchTerm)


