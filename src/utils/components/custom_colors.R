custom_colors <- function() {
  paleTurquoise <- make_style("#AFEEEE")
  aquamarine <- make_style("#7FFFD4")
  lightSteelBlue <- make_style("#B0C4DE")
  lightCoral <- make_style("#F08080")
  crimsonRed <- make_style("#DC143C")
  lightGreen <- make_style("#90EE90")
  mediumSeaGreen <- make_style("#3CB371")
  darkGreen <- make_style("#006400")
  lightSkyBlue <- make_style("#87CEFA")
  deepSkyBlue <- make_style("#00BFFF")
  dodgerBlue <- make_style("#1E90FF")
  timer <- make_style("#874949")
  darkKhaki <- make_style("#BDB76B")
  
  return(
    list(
      indicator = paleTurquoise,
      debugLight = aquamarine,
      output = darkKhaki,
      seqInit = dodgerBlue,
      funInit = deepSkyBlue,
      proInit = lightSkyBlue,
      highlight = lightSteelBlue,
      fatalError = crimsonRed,
      nonFatalError = lightCoral,
      seqSuccess = darkGreen,
      funSuccess = mediumSeaGreen,
      proSuccess = lightGreen,
      timer = timer
    )
  )
}
