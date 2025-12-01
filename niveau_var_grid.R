rm(list = ls())
source("sim_rejection_rate.R")

#### Run algorithm on 4D AR1 process, under local alternatives, including H0.

# B_func <- function(T) floor(T^(1/4))
L = 10
repetitions <- 200

set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)

seeds <- 1:200


Tlens_pwr <- c(200)
baseline_gammas <- c(0)
Bs <- c(1, 5, 10, 20, 40, 70, 100)


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



niveau_var_grid_df <- comb_rej_rate_large_obj(
  readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_1.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_5.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_20.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_1000_baseline_gamma_0_B_40.rds")
)

niveau_var_grid_df_small <- comb_rej_rate_large_obj(
  readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_1.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_5.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_20.rds"),
  readRDS("datasets/niveau/VAR/grid_Tlen_200_baseline_gamma_0_B_40.rds")
)

niveau_var_grid_df %>% 
  dplyr::select(c(alpha, rate_nonparametric, se_nonparametric, B, Tlen)) %>% 
  mutate(B = factor(B)) %>% 
  ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() +
  geom_ribbon(aes(ymin = rate_nonparametric - 1.96 * se_nonparametric,
                  ymax = rate_nonparametric + 1.96 * se_nonparametric),
              alpha = 0.2) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~B)


niveau_var_grid_df_small %>% 
  dplyr::select(c(alpha, rate_nonparametric, se_nonparametric, B, Tlen)) %>% 
  mutate(B = factor(B)) %>% 
  ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() +
  geom_ribbon(aes(ymin = rate_nonparametric - 1.96 * se_nonparametric,
                  ymax = rate_nonparametric + 1.96 * se_nonparametric),
              alpha = 0.2) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~B)


niveau_var_grid_df %>% 
  dplyr::select(c(alpha, rate_parametric_plugin, se_parametric_plugin, B, Tlen)) %>% 
  mutate(B = factor(B)) %>% 
  ggplot(aes(x = alpha, y = rate_parametric_plugin)) + 
  geom_line() +
  geom_ribbon(aes(ymin = rate_parametric_plugin - 1.96 * se_parametric_plugin,
                  ymax = rate_parametric_plugin + 1.96 * se_parametric_plugin),
              alpha = 0.2) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~B)




t1000test <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_0.rds")

t1000test$sim_rej_obj$rejection_rate_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() + 
  geom_abline(color = "darkgrey")

