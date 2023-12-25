if (!file.exists("./resources/wfo/classification.csv")) {
  tryCatch({
    if (!dir.exists("./resources/wfo")) dir.create("./resources/wfo", recursive = T)
    cat(blue("A progress window pops up here, check your taskbar \n"))
    WFO.download(save.dir = "resources/wfo", WFO.remember = TRUE, timeout = 2000)
    
  }, error = function(e) {
    cat("Error in download: ", e$message, "\n")
    cat(red("Could not download WFO Backbone. \n"))
    cat("Opening download page \n")
    
    browseURL("https://www.worldfloraonline.org/downloadData;jsessionid=D1501051E49AE20AB4B7297D021D6324")
    
  }, finally = {
    WFO_file = "./resources/wfo/classification.csv"
  }
  )
} else {
  WFO_file = "./resources/wfo/classification.csv"
}