library(patchwork)

#todo: install packages
#todo: is this needed in general?
#todo: add function descriptions to know why they do
#todo: change function name to things that makes more sense


### ------ Functions to make plots for calibration and power ------ ###
plot_rates <- function(df) {
  true_rate <- ggplot(df$rejection_rate_df, aes(x = alpha, y = rate_oracle_plugin)) +
    labs(x = "alpha", y = "Rejection rate with oracle plugin") + 
    geom_line() +
    geom_abline(color = "red") +
    theme_bw()
  
  est_rate <- ggplot(df$rejection_rate_df, aes(x = alpha, y = rate_nonparametric)) +
    labs(x = "alpha", y = "Rejection rate with nonparametrically estimated plugin") + 
    geom_line() +
    geom_abline(color = "red") +
    theme_bw()
  
  (est_rate + true_rate) +
    plot_annotation(
      title = "Estimated rejection rate",
      subtitle = "Right: based on estimated; left: based on oracle", 
      theme = theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))
    )
}

plot_scatterplot <- function(df) {
  scatter_plot <- ggplot(df$estimates) + 
    geom_point(aes(x = S_hat, y = S_oracle_plugin)) +
    geom_abline(color = "red") +
    theme_bw()
  
  scatter_plot + 
    plot_annotation(
      title = "Scatterplot of oracle vs estimated",
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
}

plot_qq_plot <- function(df) {
  qq_plot <- ggplot(df$estimates) + 
      geom_qq(aes(sample = S_hat),
              distribution = function(p) quantile(df$estimates$S_oracle_plugin, p)) +
    labs(x = "quantiles of oracle",
         y = "quantiles of estimated") +
    geom_abline(color = "red") +
    theme_bw()
  
  qq_plot + 
    plot_annotation(
      title = "QQ-plot of estimated vs oracle",
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
}

plot_histograms <- function(df) {
  histograms <- ggplot(df$estimates %>% pivot_longer(cols = c(S_hat, S_oracle_plugin))) + 
      geom_histogram(aes(x = value, y = after_stat(density)), color = "white") +
      facet_wrap(~name, labeller = as_labeller(c(S_hat = "Estimated S", S_oracle_plugin = "Oracle S"))) +
      theme_bw()
  
  histograms + 
    plot_annotation(
      title = "Histogram of test-statistics",
      subtitle = "stratified according to oracle or estimated",
      theme = theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))
    )
}

gather_plots <- function(df) {
  plot_rates(df) / (plot_scatterplot(df) + plot_qq_plot(df)) / plot_histograms(df) + 
    plot_annotation(title = "Plots for rejecetion rate",
                    theme = theme(plot.title = element_text(hjust = 0.5)))
}


### ------ functions to make plots to assess remainder terms ------ ###

plot_remainder_histograms <- function(df) {
  histograms <- ggplot(df$remainders %>% pivot_longer(cols = c(R1, R2, R3))) + 
    geom_histogram(aes(x = value, y = after_stat(density)), color = "white") + 
    facet_wrap(~name, labeller = as_labeller(
      c(R1 = "R1", R2 = "R2", R3 = "R3"))) + 
    theme_bw()
  
  histograms + 
    plot_annotation(
      title = "Distribution of remainders",
      subtitle = paste0("Length of chain: ", df$Tlen),
      theme = theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))
    )
}

#not sure this works as intended... difficult to analyze with small sample size
plot_remainder_decays <- function(...) {
  arguments <- list(...)
  number_dfs <- length(arguments)
  
  dfs <- vector("list", length = number_dfs)

  for( i in 1:number_dfs) {
    arg_i <- arguments[[i]]
    dfs[[i]] <- arg_i$remainder %>% mutate(Tlen = arg_i$Tlen)
  }
  
  combined_data_frame <- bind_rows(dfs)
  # browser()
  
  summarized_df <- combined_data_frame %>%
    group_by(Tlen) %>%
    mutate(
      R1_rate = R1 * sqrt(log(Tlen)),
      R2_rate = R2 * sqrt(log(Tlen)),
      R3_rate = R3 * sqrt(log(Tlen))
    ) %>%
    summarise(
      mean_R1_rate = mean(R1_rate, na.rm = TRUE),
      mean_R2_rate = mean(R2_rate, na.rm = TRUE),
      mean_R3_rate = mean(R3_rate, na.rm = TRUE)
    ) %>%
    ungroup()
  
  summarized_df
}



