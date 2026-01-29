### Code used for MSc presentation ###

# source relevant files #
source("simulate_AR_process.R")
source("estimate_test.R")
source("sim_crit_value.R")
source("test_hypothesis.R")
source("simulate_parameters.R")
source("helper_functions.R")
source("conditional_distributions.R")

#true conditional characteristic functions
remainder_true_ccfs = list(
  true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
  true_psi = function(x, u, A, t) char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
)


library(patchwork)
library(tidyverse)

#create A matrix

#different e-functions to use
e_flat <- function(d) rep(1/sqrt(d), d)
e_sparse <- function(d) c(1, rep(0, d - 1))
e_random <- function(d) rnorm(d)

A <- matrix(c(0.21, 0.94, -0.65, -0.05, 
              0.45, 0.75, -0.24, 0.28, 
              0.31, 0.01, -0.39, 0.43, 
              -0.42, -0.61, -0.35, 0.68), 
            4, 4) 
A

# set.seed(420)
# A <- runif(16, -1, 1) %>% matrix(4, 4) %>% round(2)

# A <- AR1_matrix(d = 4, gamma = 0)
A_DBMLI <- A
A_DBMLI[1, -1] <- 0 #DBMLI satisfied
(A_DBMLI %>% eigen())$values %>% abs() %>% max() 

#datalength, # of evaluation points and # of cross-fit-folds
Tlen <- 1000
L <- 10

B <- ceiling(Tlen^(49/100))

# evaluation points
mu <- matrix(rnorm(B), ncol = 1)
nu <- matrix(rnorm(B * (ncol(A) - 1)), ncol = (ncol(A) - 1))

#generate data 

#don't fuck with this :)
data_DBMLI <- simulate_AR_process(Tlen, A_DBMLI) #DBMLI satisfied
data <- simulate_AR_process(Tlen, A) #DBMLI not satisfied

# saveRDS(data, file = "datasets/test/data1000-6.rds")
# saveRDS(data_DBMLI, file = "datasets/test/data_DBMLI1000-6.rds")
# saveRDS(mu, file = "datasets/test/mu6.rds")
# saveRDS(nu, file = "datasets/test/nu6.rds")


#make time series plots
p_ts_DBMLI <- make_ts_plot(data_DBMLI, title = "data 1")
p_ts <- make_ts_plot(data, title = "data 2")

time_series_plot <- p_ts_DBMLI + p_ts
time_series_plot

test_hypothesis(data_DBMLI, L, B, mu, nu, 0.05)
test_hypothesis(data, L, B, mu, nu, 0.05)

#make residual plots

#residuals
residuals <- data.frame( #DBMLI not satisfied
  residuals = data[, "X"][-1] - A[1,1] * data[-nrow(data), "X"],
  Y1        = data[-nrow(data), "Y1"]
)

residuals_DBMLI <- data.frame( #DBMLI satisfied
  residuals = data_DBMLI[, "X"][-1] - A_DBMLI[1,1] * data_DBMLI[-nrow(data_DBMLI), "X"],
  Y1        = data_DBMLI[-nrow(data_DBMLI), "Y1"]
)

#residual plots
p_res <- ggplot(residuals) +  #DBMLI not satisfied
  geom_point(aes(x = Y1, y = residuals), size = 0.75) + 
  labs(
    x = expression(Y[t]^(1)),
    y = expression(R[t+1])
  ) +
  theme_bw()

p_res_DBMLI <- ggplot(residuals_DBMLI) + 
  geom_point(aes(x = Y1, y = residuals), size = 0.75) + 
  labs(
    x = expression(Y[t]^(1)),
    y = expression(R[t+1])
  ) +
  theme_bw()

#combined ts plot and residual plot
time_series_plot <- p_ts_DBMLI + p_ts

# ggsave("images/presentation/time_series_plots.eps",
#        width = 500,
#        height = 320,
#        units = "px",
#        scale = 4,
#        plot = time_series_plot)

test_hypothesis(data_DBMLI, L, B, mu, nu, 0.05)
test_hypothesis(data, L, B, mu, nu, 0.05)

time_series_plot


combined_plot <- (p_ts_DBMLI | p_ts) /
  (p_res_DBMLI | p_res) 
combined_plot

# ggsave("images/presentation/residuals.eps",
#        width = 500,
#        height = 320,
#        units = "px",
#        scale = 4,
#        plot = combined_plot)


est <- estimate_stat(
  data = data, L = L, B = B,
  mu = mu, nu = nu
)

est_DBMLI <- estimate_stat(data_DBMLI, L = L, B = B, mu, nu)


test_hypothesis(data, L, B, mu, nu, 0.05) #should reject
test_hypothesis(data_DBMLI, L, B, mu, nu, 0.05) #shall not reject

dist_plot <- est$Covvar_Est %>% 
  sim_crit_draws_alt2() %>% 
  data.frame(value = .) %>% 
  ggplot() + 
  geom_histogram(aes(x = value, y = after_stat(density)), color = "white") + 
  labs(
    x = latex2exp::TeX("$\\|N(0,\\widehat{V})\\|_{\\infty}$")
  ) +
  theme_bw()
dist_plot

dist_plot_S_hat <- dist_plot +
  geom_vline(xintercept = est$S_hat,
             color = "blue", linetype = "dashed", linewidth = 0.7) +
  annotate(
    "text",
    x = est$S_hat,
    y = Inf,
    label = latex2exp::TeX(paste0("$\\widehat{S} = ", signif(est$S_hat, 4), "$")),
    vjust = 1.2,
    hjust = -0.05,
    color = "blue"
  )
dist_plot_S_hat

dist_plot_S_hat_crit_value <- dist_plot_S_hat + 
  geom_vline(xintercept = est$Covvar_Est%>% sim_crit_draws_alt2() %>% crit_from_draws(alpha = 0.05), 
             color = "red", linetype = "dashed") + 
  annotate(
    "text",
    x = est$Covvar_Est%>% sim_crit_draws_alt2() %>% crit_from_draws(alpha = 0.05) - 0.55,
    y = 1.5,
    label = latex2exp::TeX(paste0("$\\widehat{c}_{1-\\alpha} = ", signif(
      est$Covvar_Est%>% sim_crit_draws_alt2() %>% crit_from_draws(alpha = 0.05), 4), "$")),
    vjust = 1.2,
    hjust = -0.05,
    color = "red"
  )
dist_plot_S_hat_crit_value


null_dist_plot <- est_DBMLI$Covvar_Est %>% 
  sim_crit_draws_alt2() %>% 
  data.frame(value = .) %>% 
  ggplot() + 
  geom_histogram(aes(x = value, y = after_stat(density)), color = "white", bins = 40) + 
  labs(
    x = latex2exp::TeX("$\\|N(0,\\widehat{V})\\|_{\\infty}$")
  ) +
  theme_bw()
null_dist_plot

null_dist_plot_S_hat <- null_dist_plot +
  geom_vline(xintercept = est_DBMLI$S_hat,
             color = "blue", linetype = "dashed", linewidth = 0.7) +
  annotate(
    "text",
    x = est_DBMLI$S_hat - 0.55,
    y = Inf,
    label = latex2exp::TeX(paste0("$\\widehat{S} = ", signif(est_DBMLI$S_hat, 4), "$")),
    vjust = 1.2,
    hjust = -0.05,
    color = "blue"
  )
null_dist_plot_S_hat

null_dist_plot_S_hat_crit_value <- null_dist_plot_S_hat + 
  geom_vline(xintercept = est_DBMLI$Covvar_Est%>% sim_crit_draws_alt2() %>% crit_from_draws(alpha = 0.05), 
             color = "red", linetype = "dashed") + 
  annotate(
    "text",
    x = est_DBMLI$Covvar_Est%>% sim_crit_draws_alt2() %>% crit_from_draws(alpha = 0.05),
    y = Inf,
    label = latex2exp::TeX(paste0("$\\widehat{c}_{1-\\alpha} = ", signif(
      est_DBMLI$Covvar_Est%>% sim_crit_draws_alt2() %>% crit_from_draws(alpha = 0.05), 4), "$")),
    vjust = 1.2,
    hjust = -0.05,
    color = "red"
  )
null_dist_plot_S_hat_crit_value

combined_dist_plot_S_hat <- null_dist_plot_S_hat + dist_plot_S_hat
ggsave("images/presentation/dist_plot_S_hat.eps",
       width = 500, 
       height = 320,
       units = "px",
       scale = 4,
       plot = combined_dist_plot_S_hat)

combined_dist_plot_S_hat_crit_val <- null_dist_plot_S_hat_crit_value + dist_plot_S_hat_crit_value
ggsave("images/presentation/dist_plot_S_hat_critical_value.eps",
       width = 845, #500
       height = 350, #320
       units = "px",
       scale = 3,
       plot = combined_dist_plot_S_hat_crit_val)

ggsave("images/presentation/null_dist_S_hat_crit_value.eps",
       width = 500,
       height = 320,
       units = "px",
       scale = 4,
       plot = null_dist_plot_S_hat_crit_value)

#see if you can find data where DMBLI does not reject and the other does?
test <- test_hypothesis(data, L, B, mu, nu, 0.05) #should reject
test_DBMLI <- test_hypothesis(data_DBMLI, L, B, mu, nu, 0.05) #shall not reject

test$p_val %>% signif(2)
test_DBMLI$p_val %>% signif(2)



