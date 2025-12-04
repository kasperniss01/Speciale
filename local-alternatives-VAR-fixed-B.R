### local alternatives for VAR ### 

## Y is 3D - fixed B ##

source("sim_rejection_rate.R")
source("helper_functions.R")

#structure
L <- 10
repetitions <- 200
set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)

#seeds
seeds <- 1:200

# niveau analyse hvor cross-fitting undersøges
# fast B = 10, varierende Tlen

Bs <- c(10)
Tlens_pwr <- c(100, 500, 1000, 2000)
baseline_gammas <- c(0, 3, 5, 10, 15)

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    for (B in Bs) {
      
      if(Tlen == 100 & baseline_gamma > 5) next
      
      if (baseline_gamma %in% c(0, 5) & Tlen != 2000) next
      
      print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
      
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
        repetitions = repetitions
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
              file = paste0("datasets/new_sims/VAR/_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B, ".rds"))
    
    }
  }
}


## plotting results
local_alternative_VAR_fixed_B_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/_Tlen_100_baseline_gamma_3_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_3_B_10.rds"),
  # readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_5_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_10_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_15_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_3_B_10.rds"),
  # readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_5_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_10_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_15_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_3_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_5_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_10_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_15_B_10.rds")
)

#full plot for local alternative
local_alternative_VAR_fixed_B_full_plot <- local_alternative_VAR_fixed_B_df %>% 
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
local_alternative_VAR_fixed_B_full_plot

local_alternative_VAR_fixed_B_alpha_0.05_plot <- local_alternative_VAR_fixed_B_df %>% 
  filter(alpha <= 0.0475 & alpha <= 0.0525) %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin), 
               values_to = "rate",
               names_to = "method") %>% 
  ggplot(aes(x = Tlen, y = rate, color = method)) + 
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
local_alternative_VAR_fixed_B_alpha_0.05_plot

local_alt_fast_B_alpha_0.05_plot <- local_alternative_VAR_fixed_B_df %>% filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
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
