### CIR local alternative 1 ###

source("sim_rejection_rate.R")
source("simulate_parameters.R")

L <- 10
repetitions <- 200

#Y is 3-dimensional
d <- 3
e <- 1/sqrt(d) * rep(1, d)

#seeds for param generation and datasets
seed <- 420
seeds <- 1:200

# Tlens <- c(seq(100, 1000, 100), 1000, 2000, 3000, 4000, 5000)
Tlens <- c(2000, 3000, 4000, 5000) #har kørt til og med 2000 og alle gammaer - mangler resten af T
baseline_gammas <- c(1, 10, 15, 20, 30)
Bs <- c(10)

#initial theta values
theta_base <- CIR_param(d = d, seed = seed, gamma = 0) 

for(Tlen in Tlens){
  for(baseline_gamma in baseline_gammas){
    for(B in Bs) {
      
      print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma, ", B = ", B))
      
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
        B = B,
        DGP = "CIR", 
        parameters = list(
          theta = theta_local
        ),
        alphas = seq(0.005, 1, 0.005),
        repetitions = repetitions
      )
      
      
      sim_temp <- list(
        metadata = list(
          Tlen = Tlen,
          baseline_gamma = baseline_gamma,
          actual_gamma = gamma,
          theta = theta_local, 
          # A_matrix = Alocal,
          B = B,
          L = L,
          repetitions = repetitions
        ),
        sim_rej_obj = simulate_temp
      )
      
      saveRDS(sim_temp, file = paste0("datasets/niveau/CIR/grid_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B, ".rds"))
      
    }
  }
}