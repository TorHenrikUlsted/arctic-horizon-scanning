

// Rename this function to something that makes more sense
// (I lack the domain knowledge)
//
// Add error handling if you want by checking data.status
// Should be = 200 if the request was OK 
const search = async (group, searchTerm) => {
  const response = await fetch(`https://api.gbif.org/v1/${group}/search?${searchTerm}`)
  const data = await response.json()

  console.log(data)
}

search('occurrence', 'year=1800,1899')

