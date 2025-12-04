### Cross-fitting - niveau - VAR - fixed B ###
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

Tlens_pwr <- c(30, 100, 500, 1000)
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
              file = paste0("datasets/new_sims/VAR/L_", L, "_Tlen_", Tlen, "_B_", B, ".rds"))
      
    }
  }
}

#færdig med ingen cross-fit
L <- 10

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
              file = paste0("datasets/new_sims/VAR/L_", L, "_Tlen_", Tlen, "_B_", B, ".rds"))
      
    }
  }
}


# plotting the result - no cross-fitting
# no_cross_fit_30 <- readRDS("datasets/niveau/VAR/single_crossfit_Tlen_200_baseline_gamma_0_B_10.rds")
# no_cross_fit_200 <- readRDS("datasets/niveau/VAR/single_crossfit_Tlen_200_baseline_gamma_0_B_10.rds")
# no_cross_fit_500 <- readRDS("datasets/niveau/VAR/single_crossfit_Tlen_500_baseline_gamma_0_B_10.rds")
# no_cross_fit_2000 <- readRDS("datasets/niveau/VAR/single_crossfit_Tlen_2000_baseline_gamma_0_B_10.rds")


no_cross_fit_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/L_1_Tlen_30_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_1_Tlen_100_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_1_Tlen_500_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_1_Tlen_1000_B_10.rds")
  )
cross_fit_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/L_10_Tlen_30_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_100_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_500_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_B_10.rds")
)

assess_cross_fit_df <- dplyr::bind_rows(
  no_cross_fit_df %>% dplyr::mutate(crossfit = "Without cross-fitting"),
  cross_fit_df %>% dplyr::mutate(crossfit = "With cross-fitting")
)

asses_cross_fit_plot <- assess_cross_fit_df %>% 
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

asses_cross_fit_plot
