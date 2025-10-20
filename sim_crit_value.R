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

sim_crit_draws <- function(covvar_est, nrep = 1e5) {
  m <- nrow(covvar_est)
  
  R <- t(chol(covvar_est))                                
  Z <- matrix(rnorm(m * nrep), m, nrep)
  V <- R %*% Z
  
  col_max_abs <- function(M) apply(M, 2, function(v) max(abs(v)))
  col_max_abs(V)
}

crit_from_draws <- function(draws, alpha) {
  q <- 1 - alpha
  unname(quantile(draws, q))
}