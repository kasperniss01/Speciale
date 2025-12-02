### R functionality to simulate critical value for S_hat

sim_crit_draws <- function(covvar, nrep = 1e5) {
  # covvar is a covvariance matrix
  # nrep is number of replications for quantile determination
  
  #returns vector of length(nrep) of samples from ||N(0, covvar)||_infty
  
  m <- nrow(covvar)
  
  R <- t(chol(covvar))                                
  Z <- matrix(rnorm(m * nrep), m, nrep)
  V <- R %*% Z

  Rfast::colMaxs(abs(V), value = TRUE)
}

sim_crit_draws_alt <- function(all_gamma, nrep = 1e5) {
  # all_gamma is (T - L) x (2 * B) one row pr. gamma, on col pr. B/Re/Im,{ (T-L) = n }
  # nrep is number of replications for quantile determination
  
  #returns vector of length(nrep) of samples from ||N(0, covvar)||_infty
  
  n <- nrow(all_gamma) 
  
  Zalt <- matrix(rnorm(n * nrep), nrep, n)/sqrt(n)
  
  Valt <- Zalt %*% all_gamma
  
  Rfast::rowMaxs(abs(Valt), value = TRUE)
}

crit_from_draws <- function(draws, alpha) {
  #draws a vector of draws from a distribution
  #alpha is a vector of quantiles
  
  #returns a vector of critical values from a distribution
  
  q <- 1 - alpha
  unname(quantile(draws, q))
}

# deprecated?
# sim_crit_value <- function(n = 1e3, covvar_est, alpha = 0.05) {
#   B <- nrow(covvar_est) / 2
#   dist <- replicate(n, {
#     Z <- rnorm(2 * B)
#     
#     root_covvar_est <- t(chol(covvar_est))
#     max(abs(root_covvar_est %*% Z))
#   })
#   
#   unname(quantile(dist, 1 - alpha/2))
# }
