library(patchwork)

plot_rates <- function(df) {
  true_rate <- ggplot(df$rejection_rate_df, aes(x = alpha, y = rate_true)) +
    labs(x = "alpha", y = "Estimated (true) rejection rate") + 
    geom_line() +
    geom_abline(color = "red") +
    theme_bw()
  
  est_rate <- ggplot(df$rejection_rate_df, aes(x = alpha, y = rate)) +
    labs(x = "alpha", y = "Estimated rejection rate") + 
    geom_line() +
    geom_abline(color = "red") +
    theme_bw()
  
  (est_rate + true_rate) +
    plot_annotation(
      title = "Estimated rejection rate",
      subtitle = "Left: based on oracle; right: based on estimated", 
      theme = theme(plot.title = element_text(hjust = 0.5),
                    plot.subtitle = element_text(hjust = 0.5))
    )
}

plot_scatterplot <- function(df) {
  scatter_plot <- ggplot(df$estimates) + 
    geom_point(aes(x = S_hat, y = S_true)) +
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
              distribution = function(p) quantile(df$estimates$S_true, p)) +
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
  histograms <- ggplot(df$estimates %>% pivot_longer(cols = c(S_hat, S_true))) + 
      geom_histogram(aes(x = value, y = after_stat(density)), color = "white") +
      facet_wrap(~name, labeller = as_labeller(c(S_hat = "Estimated S", S_true = "Oracle S"))) +
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


