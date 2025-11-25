rm(list = ls())
source("sim_rejection_rate.R")
source("simulate_parameters.R")

#to initialize .Random.seed
B_func <- function(T) floor(T^(1/4))
L <- 10
repetitions <- 200
seed <- 420
d <- 3 #dimension of Y is 3
e <- 1/sqrt(d) * rep(1, d)

# Tlens_pwr <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
#                2000, 3000, 4000, 5000)

Tlens_pwr <- c(1000, 2000, 3000, 4000, 5000)
baseline_gammas <- c(1, 3, 5, 7, 9, 11)

#initial theta values
theta_base <- CIR_param(d = d, seed = seed, gamma = 0) 


for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    
    print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
    
    gamma <- baseline_gamma / sqrt(Tlen)
    
    theta_base <- theta_base
    theta2_local <- theta_base$theta2
    theta2_local[1, -1] <- gamma * e
    
    theta_local <- list(theta1 = theta_base$theta1,
                        theta2 = theta2_local,
                        theta3 = theta_base$theta3)
    # Alocal <- A
    # Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B_func(Tlen),
      DGP = "CIR", 
      parameters = list(
        theta = theta_local
      ),
      alphas = seq(0.005, 1, 0.005),
      # remainder_true_ccfs = list(
      #   true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
      #   true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
      # ),
      # parametric = TRUE,
      repetitions = repetitions
    )
    
    
    sim_temp <- list(
      metadata = list(
        Tlen = Tlen,
        baseline_gamma = baseline_gamma,
        actual_gamma = gamma,
        theta = theta, 
        # A_matrix = Alocal,
        B = B_func(Tlen),
        L = L,
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    saveRDS(sim_temp, file = paste0("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, ".rds"))
    
    
  }
}
