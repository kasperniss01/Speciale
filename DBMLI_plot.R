library(patchwork)

set.seed(420)
A <- matrix(runif(16, -1, 1), ncol = 4) %>% round(2)

Tlen <- 100

Ddata <- simulate_AR_process(Tlen, A)
Dplot <- plot(Ddata, main = "")

Dresiduals <- data.frame(
  residuals = Ddata[, "X"][-1] - A[1,1] * Ddata[-nrow(Ddata), "X"],
  Y1 = Ddata[-nrow(Ddata), "Y1"]
)

ggplot(Dresiduals) + 
  geom_point(aes(x = Y1, y = residuals)) + 
  labs(
    x = expression(Y^1),
    y = expression(R[t+1])
  ) + 
  theme_bw()


A[1, -1] <- 0
Idata <- simulate_AR_process(Tlen, A)
Iplot <- plot(Idata, main = "")

Iresiduals <- data.frame(
  residuals = Idata[, "X"][-1] - A[1,1] * Idata[-nrow(Idata), "X"],
  Y1 = Idata[-nrow(Idata), "Y1"]
)

ggplot(Iresiduals) + 
  geom_point(aes(x = Y1, y = residuals)) + 
  labs(
    x = expression(Y^1),
    y = expression(R[t+1])
  ) + 
  theme_bw()




###
Tlen <- 100
set.seed(420)
A <- matrix(runif(4, -1, 1), ncol = 2) %>% round(2)

Ddata <- simulate_AR_process(Tlen, A)
Dresiduals <- data.frame(
  residuals = Ddata[, "X"][-1] - A[1,1] * Ddata[-nrow(Ddata), "X"],
  Y1        = Ddata[-nrow(Ddata), "Y1"]
)

A[1, -1] <- 0
Idata <- simulate_AR_process(Tlen, A)

Iresiduals <- data.frame(
  residuals = Idata[, "X"][-1] - A[1,1] * Idata[-nrow(Idata), "X"],
  Y1        = Idata[-nrow(Idata), "Y1"]
)

## --- helper: ggplot version of the multivariate time series --------
make_ts_plot <- function(dat, title) {
  df <- as.data.frame(dat)
  
  df$Time <- seq_len(nrow(df)) / frequency(dat)

  
  long <- pivot_longer(df, cols = -Time,
                       names_to = "series", values_to = "value")
  
  long$series <- factor(
    long$series,
    levels = c("X", "Y1", "Y2", "Y3"),
    labels = c("X", "Y^(1)", "Y^(2)", "Y^(3)")
  )
  
  ggplot(long, aes(x = Time, y = value)) +
    geom_line(linewidth = 0.25) +
    facet_grid(
      series ~ ., 
      # scales = "free_y", 
      switch = "y",
      labeller = label_parsed       # <- parse "Y^1" etc as expressions
    ) +
    labs(x = "Time", y = NULL, title = title) +
    theme_bw() +
    theme(
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0),
      plot.title = element_text(hjust = 0.5),
      panel.spacing.y = unit(0.15, "lines")
    )
}

p_ts_I <- make_ts_plot(Idata, "DBMLI true")
p_ts_D <- make_ts_plot(Ddata, "DBMLI violated")

p_ts_D <- p_ts_D +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

p_res_D <- p_res_D +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank())

## --- residual plots ------------------------------------------------
p_res_I <- ggplot(Iresiduals) + 
  geom_point(aes(x = Y1, y = residuals)) + 
  labs(
    x = expression(Y[t]^(1)),
    y = expression(R[t+1])
  ) +
  theme_bw()

p_res_D <- ggplot(Dresiduals) + 
  geom_point(aes(x = Y1, y = residuals)) + 
  labs(
    x = expression(Y[t]^(1)),
    y = expression(R[t+1])
  ) +
  theme_bw()

## --- 2x2 layout: top = time series, bottom = residuals -------------
combined_plot <- (p_ts_I | p_ts_D) /
  (p_res_I | p_res_D) 
combined_plot

ggsave("images/BMLI_violations.eps", 
       width = 500,
       height = 500,
       units = "px",
       scale = 4,
       plot = combined_plot)



# gridExtra::grid.arrange(list(Dplot, resDplot, Iplot, resIplot), nrow = 2, ncol = 2)
theta_base <- CIR_param(d = 3, seed = 420, gamma = 0) 
drift_base <- make_CIR_drift(theta_base$theta1, theta_base$theta2)
diffusion_base <- make_CIR_diffusion(theta_base$theta3)
SDE_data <- simulate_sde(100, drift = drift_base, diffusion = diffusion_base, Z0 = rep(1, 4))


full_path <- SDE_data$whole_path %>% make_ts_plot(title = "Full path") + ylim(0, 2.5)
disc_path <- SDE_data$discretized_path %>% make_ts_plot(title = "Discretized path") + ylim(0, 2.5)

combined_paths <- full_path | disc_path
combined_paths
ggsave("images/full_and_discretized_path.eps", 
       width = 500,
       height = 300,
       units = "px",
       scale = 4,
       plot = combined_paths)
