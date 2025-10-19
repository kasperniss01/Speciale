### document for simulating rejection rate 
### right now only works for 2D-AR(1) processes!!!

source("simulate_AR_process.R")
source("estimate_test.R")
source("sim_crit_value.R")

sim_rej_rate <- function(Tlen, n, L, B,
                         A_matrix, 
                         alpha,
                         repetitions = 500,
                         verbose = TRUE) {
  
  mu <- matrix(rnorm(B), ncol = 1)
  nu <- matrix(rnorm(B), ncol = 1)   # d = 1 here
  
  reject <- logical(repetitions)
  
  for (i in seq_len(repetitions)) {
    data <- simulate_2D_AR1_process(Tlen, A_matrix, verbose = FALSE)
    
    est <- estimate_stat(
      data = data, n = n, L = L, B = B,
      mu = mu, nu = nu
    )
    
    S_hat <- est$S_hat
    crit  <- sim_crit_value(covvar_est = est$Covvar_Est, alpha = alpha)  
    
    reject[i] <- (S_hat > crit)
    if (verbose && i %% 10 == 0) message(i)
  }
  
  mean(reject)
}


Tlen <- 1000
n <- 200
L <- 5
B <- 10

A <- matrix(c(-0.4, 0, -0.3, 0.8), nrow = 2, byrow = T)

sim_rej_rate(Tlen, n, L, B, A, 0.05)
