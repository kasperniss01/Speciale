source("sim_rejection_rate.R")

L <- 10
repetitions <- 200

set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)

seeds <- 1:200

Tlens_pwr <- c(200, 500)
baseline_gammas <- c(1, 5, 10, 15)
Bs <- c(10)

# cat("Lokale alternativer, fast B =", B, "varierende T", "forskellige gamma_0 = ", baseline_gammas)

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