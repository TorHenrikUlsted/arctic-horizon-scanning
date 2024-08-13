# Create a list to store the timers
timers <- list()

start_timer <- function(id) {
  # If a timer with the given ID already exists, add an interval at the end of the ID
  interval <- 1
  while(!is.null(timers[[paste0(id, " ", interval)]])) {
    interval <- interval + 1
  }
  id <- paste0(id, "_", interval)
  
  timers[[id]] <<- list(start_time = Sys.time())
  vebcat("Time tracking started for", id, color = "timer")
  
  return(id)
}

calculate_time <- function(time) {
  days <- floor(time / (24*60*60))
  time <- time - (days*24*60*60)
  
  hours <- floor(time / (60*60))
  time <- time - (hours*60*60)
  
  minutes <- floor(time / 60)
  seconds <- time - (minutes*60)
  
  return(list(
    days = days,
    hrs = hours,
    mins = minutes,
    secs = seconds
  ))
}

end_timer <- function(id) {
  if(is.null(timers[[id]])) {
    vebcat("Error: Timer", id, "was not started.", color = "timer")
  } else {
    end_time <- Sys.time()
    time_diff <- as.numeric(difftime(end_time, timers[[id]]$start_time, units = "secs"))
    
    time <- calculate_time(time_diff)
    
    vebcat("Time elapsed for timer", id, ": ", time$days, "day(s)", time$hrs, "hour(s)", time$mins, "minute(s)", round(time$secs, 2), "second(s)", color = "timer")
    
    dt <- data.table(
      id = id, 
      days = time$days, 
      hours = time$hrs, 
      minutes = time$mins, 
      seconds = round(time$secs, 2))
    
    return(time_diff)
  }
}

calculate_etc <- function(timer.res, cores = 1, data.length = 1) {
  etc <- (timer.res / cores) * data.length
  
  time <- calculate_time(etc)
  
  vebcat("Estimated time to completion:", time$days, "day(s)", time$hrs, "hour(s)", time$mins, "minute(s)", round(time$secs, 2), "second(s)", color = "timer")
  
  return(time)
}

print_all_timers <- function() {
  # Initialize an empty data frame to store the timer information
  timer_df <- data.frame(id = character(), days = numeric(), hours = numeric(), minutes = numeric(), seconds = numeric())
  
  for(id in names(timers)) {
    timer_info <- end_timer(id)
    if(!is.null(timer_info)) {
      timer_df <- rbind(timer_df, timer_info)
    }
  }
  
  # Write the data frame to a CSV file
  if (!dir.exists("./outputs/post_process/")) dir.create("./outputs/post_process/", recursive = T)
  fwrite(timer_df, "./outputs/post_process/timers.csv", row.names = FALSE, bom = TRUE)
}

