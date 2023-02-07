var XMLHttpRequest = require('xhr2');
var request = new XMLHttpRequest();
var group = "/occurence";
var searchKey = "/Pan troglodytes";

request.open('GET', 'https://api.gbif.org/v1' + {group} + {searchKey}, true)
request.onload = function () {
  // Begin accessing JSON data here
  var data = JSON.parse(this.response)

  if (request.status >= 200 && request.status < 400) {
    data.forEach((data) => {
      console.log(data)
    })
  } else {
    console.log('error')
  }
}

request.send();

