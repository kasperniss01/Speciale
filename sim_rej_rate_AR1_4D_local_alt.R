source("sim_rejection_rate.R")

#### Run algorithm on 4D AR1 process, under local alternatives, including H0.

Tlens_pwr <- c(1000, 2000, 3000, 4000, 5000)
baseline_gammas <- c(0,1,3,5,7)
B_func <- function(T) floor(T^(1/4))
L = 10
repetitions <- 200

set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
    
    gamma = baseline_gamma / sqrt(Tlen)
    
    Alocal <- A
    Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B_func(Tlen),
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
        B = B_func(Tlen),
        L = L,
        repetitions = reppetitions
      ),
      sim_rej_obj = simulate_temp,
    )
    
    saveRDS(sim_temp, file = paste0("datasets/sim_rej_rate_AR1_4D_local_alt_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, ".rds"))
    
 
  }
}



