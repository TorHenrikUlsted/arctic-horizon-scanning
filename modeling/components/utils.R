# Libraries
library(dplyr)
## Get WorldFlora Library and load classification material
library(WorldFlora)
## Download and remember data
if (!file.exists("resources/classification.csv")) {
  WFO.download(save.dir = "resources", WFO.remember = TRUE)
} else {
  WFO_file = "resources/classification.csv"
}

## Time tracker
format_elapsed_time = function(start_time, end_time) {
  ### Calculate elapsed time in hours and minutes
  elapsed_time_hours = as.numeric(difftime(end_time, start_time, units = "hours"))
  elapsed_time_minutes = as.numeric(difftime(end_time, start_time, units = "mins"))
  
  ### Format elapsed time
  formatted_elapsed_time = paste(round(elapsed_time_hours, 2), "hours (", round(elapsed_time_minutes, 2), "minutes)")
  
  return(formatted_elapsed_time)
}