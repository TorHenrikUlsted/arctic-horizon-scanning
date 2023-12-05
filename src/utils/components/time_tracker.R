# Create a list to store the timers
timers <- list()

start_timer <- function(id) {
  # If a timer with the given ID already exists, add an interval at the end of the ID
  interval <- 1
  while(!is.null(timers[[paste0(id, "_", interval)]])) {
    interval <- interval + 1
  }
  id <- paste0(id, "_", interval)
  
  timers[[id]] <<- list(start_time = Sys.time())
  cat(magenta(paste("Time tracking started for timer", id, "\n")))
  
  return(id)
}

end_timer <- function(id) {
  if(is.null(timers[[id]])) {
    cat(magenta(paste("Error: Timer", id, "was not started.\n")))
  } else {
    end_time <- Sys.time()
    time_diff <- as.numeric(difftime(end_time, timers[[id]]$start_time, units = "secs"))
    
    days <- floor(time_diff / (24*60*60))
    time_diff <- time_diff - (days*24*60*60)
    
    hours <- floor(time_diff / (60*60))
    time_diff <- time_diff - (hours*60*60)
    
    minutes <- floor(time_diff / 60)
    seconds <- time_diff - (minutes*60)
    
    cat(magenta(paste("Time elapsed for timer", id, ": ", days, "days", hours, "hours", minutes, "minutes", round(seconds, 2), "seconds\n")))
    
    return(data.frame(id = id, days = days, hours = hours, minutes = minutes, seconds = round(seconds, 2)))
  }
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