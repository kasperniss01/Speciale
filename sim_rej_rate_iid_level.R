rm(list = ls())

source("sim_rejection_rate.R")

simulate_iid_gaussian_ts <- function(Tlen) {

  d <- 3
  
  mu <- rep(0, d+1)
  Sigma <- diag(d+1)
  
  data <- MASS::mvrnorm(Tlen, mu = mu, Sigma = Sigma)
  
  colnames(data) <- c("X", paste0("Y", 1:d))

  
  return(data %>% ts)  
}




Tlens_iid <- c(100)
for(Tlen in Tlens_iid){

    print(paste0("Tlen = ", Tlen))
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B_func(Tlen),
      parameters = list(
        "distribution" = simulate_iid_gaussian_ts
      ),
      DGP = "IID",
      alphas = seq(0.005, 1, 0.005),
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
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    saveRDS(sim_temp, file = paste0("datasets/iid_level_analysis/sim_rej_rate_IID_Gaussian_4D_Tlen_", Tlen, ".rds"))
    
    
}

