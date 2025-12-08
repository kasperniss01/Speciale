### CIR increasing B - level and local alternatives ### 

source("sim_rejection_rate.R")
source("simulate_parameters.R")

#setup
L <- 10
repetitions <- 200
B_func <- function(T) ceiling(T^(49 / 100))

#Y is 3-dimensional
d <- 3
e <- 1/sqrt(d) * rep(1, d)

#seeds for param generation and datasets
seed <- 420
seeds <- 1:200

# Tlens <- c(seq(100, 1000, 100), 1000, 2000, 3000, 4000, 5000)
Tlens <- c(100, 500, 1000, 2000)
baseline_gammas <- c(0, 1, 10, 20, 30)

#initial theta values
theta_base <- CIR_param(d = d, seed = seed, gamma = 0) 

for(Tlen in Tlens){
  for(baseline_gamma in baseline_gammas){
      
      print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
      
      gamma <- baseline_gamma / sqrt(Tlen)
      
      theta_base <- theta_base
      theta2_local <- theta_base$theta2
      theta2_local[1, -1] <- gamma * e
      
      theta_local <- list(theta1 = theta_base$theta1,
                          theta2 = theta2_local,
                          theta3 = theta_base$theta3)
      
      simulate_temp <- sim_rej_rate(
        Tlen = Tlen,
        L = L,
        B = B_func(Tlen),
        DGP = "CIR", 
        parameters = list(
          theta = theta_local
        ),
        alphas = seq(0.005, 1, 0.005),
        repetitions = repetitions
        # random_seeds = seeds
      )
      
      
      sim_temp <- list(
        metadata = list(
          Tlen = Tlen,
          baseline_gamma = baseline_gamma,
          actual_gamma = gamma,
          theta = theta_local, 
          # A_matrix = Alocal,
          B = B_func(Tlen),
          L = L,
          repetitions = repetitions
        ),
        sim_rej_obj = simulate_temp
      )
      
      saveRDS(sim_temp, 
              file = paste0("datasets/new_sims/CIR/_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B_func(Tlen), ".rds"))
      
    
  }
}
