rm(list = ls())
source("sim_rejection_rate.R")
library(patchwork)

L <- 10
repetitions <- 200

set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)

seeds <- 1:200

# niveau analyse hvor cross-fitting undersøges
# fast B = 10, varierende Tlen
Tlens_pwr <- c(200, 500, 2000)
baseline_gammas <- c(0)
Bs <- c(10)
L <- 1

cat("Cross-fitting, L =", L)

#running analysis
for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    for (B in Bs) {
      
      print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma, ", B = ", B))
      
      gamma = baseline_gamma / sqrt(Tlen)
      
      Alocal <- A
      Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
      
      simulate_temp <- sim_rej_rate(
        Tlen = Tlen,
        L = L,
        B = B,
        parameters = list(
          A_matrix = Alocal
        ),
        alphas = seq(0.005, 1, 0.005),
        remainder_true_ccfs = list(
          true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
          true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
        ),
        parametric = TRUE,
        repetitions = repetitions,
        random_seeds = seeds
      )
      
      
      sim_temp <- list(
        metadata = list(
          Tlen = Tlen,
          baseline_gamma = baseline_gamma,
          actual_gamma = gamma,
          A_matrix = Alocal,
          B = B,
          L = L,
          repetitions = repetitions
        ),
        sim_rej_obj = simulate_temp
      )
      
      saveRDS(sim_temp, 
              file = paste0("datasets/niveau/VAR/single_crossfit_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B, ".rds"))
      
    }
  }
}

# plotting the result - no cross-fitting
no_cross_fit_200 <- readRDS("datasets/niveau/VAR/single_crossfit_Tlen_200_baseline_gamma_0_B_10.rds")
no_cross_fit_500 <- readRDS("datasets/niveau/VAR/single_crossfit_Tlen_500_baseline_gamma_0_B_10.rds")
no_cross_fit_2000 <- readRDS("datasets/niveau/VAR/single_crossfit_Tlen_2000_baseline_gamma_0_B_10.rds")

no_cross_fit_df <- comb_rej_rate_large_obj(no_cross_fit_200,
                                        no_cross_fit_500,
                                        no_cross_fit_2000)

#combining L = 1 plots
without_cross_fit_plot <- no_cross_fit_df %>% 
  dplyr::select(-starts_with("se"), -rate_parametric) %>% 
  pivot_longer(cols = c(rate_nonparametric),
               values_to = "rate",
               names_to = "method",
               names_prefix = "rate_") %>% 
  ggplot(aes(x = alpha, y = rate)) +
  geom_line() +
  labs(
    title = "Without cross-fitting",
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~Tlen + B + baseline_gamma,
             ncol = 1,
             labeller = label_bquote(
               Tlen == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  coord_fixed() + 
  theme_bw() + 
  theme(
    strip.background = element_rect(
      fill  = "white",
      colour = NA),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5))
without_cross_fit_plot #actually looks decent??

#færdig med cross-fit analyse
L <- 10 #sættes tilbage til det oprindelige

# niveau analyse
# fast B = 10, varierende Tlen
Tlens_pwr <- c(2000, 200, 500, 1000, 3000, 4000, 5000)
baseline_gammas <- c(0)
Bs <- c(10)

cat("Niveau-analyse, fast B =", B, "varierende T")

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    for (B in Bs) {
    
    print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma, ", B = ", B))
    
    gamma = baseline_gamma / sqrt(Tlen)
    
    Alocal <- A
    Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B,
      parameters = list(
        A_matrix = Alocal
      ),
      alphas = seq(0.005, 1, 0.005),
      remainder_true_ccfs = list(
        true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
        true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
      ),
      parametric = TRUE,
      repetitions = repetitions,
      random_seeds = seeds
    )
    
    
    sim_temp <- list(
      metadata = list(
        Tlen = Tlen,
        baseline_gamma = baseline_gamma,
        actual_gamma = gamma,
        A_matrix = Alocal,
        B = B,
        L = L,
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    saveRDS(sim_temp, 
            file = paste0("datasets/niveau/VAR/grid_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B, ".rds"))
    
    }
  }
}

#Comparing with and without cross-fit
cross_fit_200 <- readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_10.rds")
cross_fit_500 <- readRDS("datasets/niveau/VAR/grid_Tlen_500_baseline_gamma_0_B_10.rds")
cross_fit_2000 <- readRDS("datasets/niveau/VAR/grid_Tlen_2000_baseline_gamma_0_B_10.rds")

cross_fit_df <- comb_rej_rate_large_obj(cross_fit_200,
                                           cross_fit_500,
                                           cross_fit_2000)

#combining L = 10 plots
cross_fit_plot <- cross_fit_df %>% 
  dplyr::select(-starts_with("se"), -rate_parametric) %>% 
  pivot_longer(cols = c(rate_nonparametric),
               values_to = "rate",
               names_to = "method",
               names_prefix = "rate_") %>% 
  ggplot(aes(x = alpha, y = rate)) +
  geom_line() +
  labs(
    title = "With cross-fitting",
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~Tlen + B + baseline_gamma,
             ncol = 1,
             labeller = label_bquote(
               Tlen == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  coord_fixed() + 
  theme_bw() + 
  theme(
    strip.background = element_rect(
      fill  = "white",  
      colour = NA),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)) 
cross_fit_plot 


# combined plot assessing both
combined_df <- dplyr::bind_rows(
  no_cross_fit_df  %>% dplyr::mutate(crossfit = "Without cross-fitting"),
  cross_fit_df     %>% dplyr::mutate(crossfit = "With cross-fitting")
)

asses_cross_fit <- combined_df %>% 
  dplyr::select(-starts_with("se"), -rate_parametric) %>% 
  tidyr::pivot_longer(
    cols = c(rate_nonparametric),
    values_to = "rate",
    names_to  = "method",
    names_prefix = "rate_"
  ) %>% 
  ggplot(aes(x = alpha, y = rate)) +
  geom_line() +
  geom_abline(color = "darkgrey", linetype = "dashed") +
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate"
  ) +
  facet_grid(
    crossfit ~ Tlen + B + baseline_gamma,
    labeller = label_bquote(
      cols = Tlen == .(Tlen) ~ "," ~ B == .(B)
    )
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    legend.position  = "bottom"
  )

##not a huge difference between the two, slightly better with cross-fitting

# general niveau analyse - med cross-fitting!
cross_fit_3000 <- readRDS("datasets/niveau/VAR/grid_Tlen_3000_baseline_gamma_0_B_10.rds")
cross_fit_4000 <- readRDS("datasets/niveau/VAR/grid_Tlen_4000_baseline_gamma_0_B_10.rds")
cross_fit_5000 <- readRDS("datasets/niveau/VAR/grid_Tlen_5000_baseline_gamma_0_B_10.rds")

niveau_fast_B_df <- bind_rows(cross_fit_df, 
                              comb_rej_rate_large_obj(cross_fit_3000,
                                            cross_fit_4000,
                                            cross_fit_5000))

niveau_fast_B_plot <- niveau_fast_B_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) +
  geom_line() +
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate"
  ) +
  facet_wrap(~Tlen,
             labeller = label_bquote(
               cols = Tlen == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom"
  )
niveau_fast_B_plot



# lokale alternativer
# fast B = 10, varierende Tlen, udsnit af der hvor testen holder niveau
Tlens_pwr <- c(1000, 2000, 3000, 4000, 5000)
baseline_gammas <- c(1, 5, 10, 15)
Bs <- c(10)

cat("Lokale alternativer, fast B =", B, "varierende T", "forskellige gamma_0 = ", baseline_gammas)

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    for (B in Bs) {
      
      print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma, ", B = ", B))
      
      gamma = baseline_gamma / sqrt(Tlen)
      
      Alocal <- A
      Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
      
      simulate_temp <- sim_rej_rate(
        Tlen = Tlen,
        L = L,
        B = B,
        parameters = list(
          A_matrix = Alocal
        ),
        alphas = seq(0.005, 1, 0.005),
        remainder_true_ccfs = list(
          true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
          true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
        ),
        parametric = TRUE,
        repetitions = repetitions,
        random_seeds = seeds
      )
      
      
      sim_temp <- list(
        metadata = list(
          Tlen = Tlen,
          baseline_gamma = baseline_gamma,
          actual_gamma = gamma,
          A_matrix = Alocal,
          B = B,
          L = L,
          repetitions = repetitions
        ),
        sim_rej_obj = simulate_temp
      )
      
      saveRDS(sim_temp, 
              file = paste0("datasets/niveau/VAR/grid_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B, ".rds"))
      
    }
  }
}
#T = 1000, forskellige gammaer
local_alt_1000_gamma_1 <- readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_1_B_10.rds")
local_alt_1000_gamma_5 <- readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_5_B_10.rds")
local_alt_1000_gamma_10 <- readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_10_B_10.rds")
local_alt_1000_gamma_15 <- readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_15_B_10.rds")

#T = 2000, forskellige gammaer
local_alt_2000_gamma_1 <- readRDS("datasets/niveau/VAR/grid_Tlen_2000_baseline_gamma_1_B_10.rds")
local_alt_2000_gamma_5 <- readRDS("datasets/niveau/VAR/grid_Tlen_2000_baseline_gamma_5_B_10.rds")
local_alt_2000_gamma_10 <- readRDS("datasets/niveau/VAR/grid_Tlen_2000_baseline_gamma_10_B_10.rds")
local_alt_2000_gamma_15 <- readRDS("datasets/niveau/VAR/grid_Tlen_2000_baseline_gamma_15_B_10.rds")

#T = 3000, forskellige gammaer
local_alt_3000_gamma_1 <- readRDS("datasets/niveau/VAR/grid_Tlen_3000_baseline_gamma_1_B_10.rds")
local_alt_3000_gamma_5 <- readRDS("datasets/niveau/VAR/grid_Tlen_3000_baseline_gamma_5_B_10.rds")
local_alt_3000_gamma_10 <- readRDS("datasets/niveau/VAR/grid_Tlen_3000_baseline_gamma_10_B_10.rds")
local_alt_3000_gamma_15 <- readRDS("datasets/niveau/VAR/grid_Tlen_3000_baseline_gamma_15_B_10.rds")

#T = 4000, forskellige gammaer
local_alt_4000_gamma_1 <- readRDS("datasets/niveau/VAR/grid_Tlen_4000_baseline_gamma_1_B_10.rds")
local_alt_4000_gamma_5 <- readRDS("datasets/niveau/VAR/grid_Tlen_4000_baseline_gamma_5_B_10.rds")
local_alt_4000_gamma_10 <- readRDS("datasets/niveau/VAR/grid_Tlen_4000_baseline_gamma_10_B_10.rds")
local_alt_4000_gamma_15 <- readRDS("datasets/niveau/VAR/grid_Tlen_4000_baseline_gamma_15_B_10.rds")
# 
# #T = 5000, forskellige gammaer
# local_alt_5000_gamma_1 <- readRDS("datasets/niveau/VAR/grid_Tlen_5000_baseline_gamma_1_B_10.rds")
# local_alt_5000_gamma_5 <- readRDS("datasets/niveau/VAR/grid_Tlen_5000_baseline_gamma_5_B_10.rds")
# local_alt_5000_gamma_10 <- readRDS("datasets/niveau/VAR/grid_Tlen_5000_baseline_gamma_10_B_10.rds")
# local_alt_5000_gamma_15 <- readRDS("datasets/niveau/VAR/grid_Tlen_5000_baseline_gamma_15_B_10.rds")

local_alt_df <- comb_rej_rate_large_obj(
  local_alt_1000_gamma_1,
  local_alt_1000_gamma_5,
  local_alt_1000_gamma_10,
  local_alt_1000_gamma_15,
  local_alt_2000_gamma_1,
  local_alt_2000_gamma_5,
  local_alt_2000_gamma_10,
  local_alt_2000_gamma_15,
  local_alt_3000_gamma_1,
  local_alt_3000_gamma_5,
  local_alt_3000_gamma_10,
  local_alt_3000_gamma_15,
  local_alt_4000_gamma_1,
  local_alt_4000_gamma_5,
  local_alt_4000_gamma_10,
  local_alt_4000_gamma_15
)

#full plot for local alternative
local_alt_fast_B_full_plot <- local_alt_df %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin), 
               values_to = "rate",
               names_to = "method") %>% 
  ggplot(aes(x = alpha, y = rate, color = method)) + 
  geom_line() + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") +
  ylim(0, 1) + 
  facet_grid(baseline_gamma ~ Tlen,
             labeller = label_bquote(
               rows = gamma[0] == .(baseline_gamma), 
               cols = Tlen   == .(Tlen)             
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", color = NA)
  )
local_alt_fast_B_full_plot

#local alternative plot zooming in on alpha = 5%
local_alt_fast_B_alpha_0.05_plot <- local_alt_df %>% filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
  ggplot(aes(x = Tlen, y = rate_nonparametric)) + 
  geom_line() + 
  labs(x = "Tlen",
       y = expression(paste("Estimated rejection rate (", alpha, " = 5%)"))) + 
  ylim(0, 1) +
  geom_point(size = 0.5) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~baseline_gamma,
               labeller = label_bquote(
                 cols = gamma[0] == .(baseline_gamma) ~ "," ~ B == .(B)
               )) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NA))
local_alt_fast_B_alpha_0.05_plot




#udvider lokale alternativer analysen med alle T fra niveau
Tlens_pwr <- c(200, 500)
baseline_gammas <- c(1, 5, 10, 15)
Bs <- c(10)

cat("Lokale alternativer, fast B =", B, "varierende T", "forskellige gamma_0 = ", baseline_gammas)

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    for (B in Bs) {
      
      print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma, ", B = ", B))
      
      gamma = baseline_gamma / sqrt(Tlen)
      
      Alocal <- A
      Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
      
      simulate_temp <- sim_rej_rate(
        Tlen = Tlen,
        L = L,
        B = B,
        parameters = list(
          A_matrix = Alocal
        ),
        alphas = seq(0.005, 1, 0.005),
        remainder_true_ccfs = list(
          true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
          true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
        ),
        parametric = TRUE,
        repetitions = repetitions,
        random_seeds = seeds
      )
      
      
      sim_temp <- list(
        metadata = list(
          Tlen = Tlen,
          baseline_gamma = baseline_gamma,
          actual_gamma = gamma,
          A_matrix = Alocal,
          B = B,
          L = L,
          repetitions = repetitions
        ),
        sim_rej_obj = simulate_temp
      )
      
      saveRDS(sim_temp, 
              file = paste0("datasets/niveau/VAR/grid_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B, ".rds"))
      
    }
  }
}




# 
# niveau_var_grid_df <- comb_rej_rate_large_obj(
#   readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_1.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_5.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_10.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_20.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_40.rds")
# )
# 
# niveau_var_grid_df_small <- comb_rej_rate_large_obj(
#   readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_1.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_5.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_10.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_20.rds"),
#   readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_40.rds")
# )
# 
# niveau_var_grid_df %>% 
#   dplyr::select(c(alpha, rate_nonparametric, se_nonparametric, B, Tlen)) %>% 
#   mutate(B = factor(B)) %>% 
#   ggplot(aes(x = alpha, y = rate_nonparametric)) + 
#   geom_line() +
#   geom_ribbon(aes(ymin = rate_nonparametric - 1.96 * se_nonparametric,
#                   ymax = rate_nonparametric + 1.96 * se_nonparametric),
#               alpha = 0.2) + 
#   geom_abline(color = "darkgrey", linetype = "dashed") + 
#   facet_wrap(~B)
# 
# 
# niveau_var_grid_df_small %>% 
#   dplyr::select(c(alpha, rate_nonparametric, se_nonparametric, B, Tlen)) %>% 
#   mutate(B = factor(B)) %>% 
#   ggplot(aes(x = alpha, y = rate_nonparametric)) + 
#   geom_line() +
#   geom_ribbon(aes(ymin = rate_nonparametric - 1.96 * se_nonparametric,
#                   ymax = rate_nonparametric + 1.96 * se_nonparametric),
#               alpha = 0.2) + 
#   geom_abline(color = "darkgrey", linetype = "dashed") + 
#   facet_wrap(~B)
# 
# 
# niveau_var_grid_df %>% 
#   dplyr::select(c(alpha, rate_parametric_plugin, se_parametric_plugin, B, Tlen)) %>% 
#   mutate(B = factor(B)) %>% 
#   ggplot(aes(x = alpha, y = rate_parametric_plugin)) + 
#   geom_line() +
#   geom_ribbon(aes(ymin = rate_parametric_plugin - 1.96 * se_parametric_plugin,
#                   ymax = rate_parametric_plugin + 1.96 * se_parametric_plugin),
#               alpha = 0.2) + 
#   geom_abline(color = "darkgrey", linetype = "dashed") + 
#   facet_wrap(~B)
# 
# 
# 
# 
# t1000test <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_0.rds")
# 
# t1000test$sim_rej_obj$rejection_rate_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) + 
#   geom_line() + 
#   geom_abline(color = "darkgrey")

