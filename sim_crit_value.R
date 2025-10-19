### R functionality to simulate critical value for S_hat

source("estimate_test.R")

sim_crit_value <- function(n = 1e3, covvar_est, alpha = 0.05) {
  B <- nrow(covvar_est) / 2
  dist <- replicate(n, {
    Z <- rnorm(2 * B)
    
    root_covvar_est <- t(chol(covvar_est))
    max(abs(root_covvar_est %*% Z))
  })
  
  unname(quantile(dist, 1 - alpha/2))
}