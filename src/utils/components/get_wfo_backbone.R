if (!file.exists("./resources/wfo/classification.csv")) {
  if (!file.exists("./resources/wfo/WFO_Backbone.zip")) {
    tryCatch({
      if (!dir.exists("./resources/wfo")) dir.create("./resources/wfo", recursive = T)
      cat(blue("A progress window pops up here, check your taskbar \n"))
      suppressWarnings(WFO.download(save.dir = "./resources/wfo", WFO.remember = TRUE, timeout = 2000))
    }, error = function(e) {
      cat(red("Error in download:"), e$message, "\n")
      cat(red("Could not download WFO Backbone. Manual download needed. \n"))

      cat("Opening download page \n")
      browseURL("https://www.worldfloraonline.org/downloadData;jsessionid=D1501051E49AE20AB4B7297D021D6324")
    }, finally = {
      WFO_file <- "./resources/wfo/classification.csv"
    })
  } else {
    unzip("./resources/wfo/WFO_Backbone.zip", exdir = "./resources/wfo")
    WFO_file <- "./resources/wfo/classification.csv"
  }
} else {
  WFO_file <- "./resources/wfo/classification.csv"
}
