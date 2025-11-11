### functionality to plot time series ###
plot_time_series <- function(ts_obj) {
  #ts_obj is a time_series object
  
  #returns a plot of the time series
  
  df <- as.data.frame(ts_obj)
  
  df$time <- time(ts_obj)
  
  df_long <- df %>% 
    pivot_longer(
      cols = -time,
      names_to = "Series",
      values_to = "Value"
    )
  
  ggplot(df_long, aes(x = time, y = Value, color = Series)) + 
    geom_line() + 
    labs(title = "Time series",
         x = "Time",
         y = "Value",
         color = "Series") + 
    theme_bw()
}


plot_time_series(my_X$whole_path)
plot_time_series(my_X$discretized_path)

