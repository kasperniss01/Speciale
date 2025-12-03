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

sim_crit_draws_alt2 <- function(covvar, nrep = 1e5, tol = 1e-10) {
  # all_gamma: (n x p) with p = 2B
  # returns length-nrep vector of ||N(0, covvar)||_∞ samples
  # where covvar = (1/n) t(all_gamma) %*% all_gamma
  
  #n <- nrow(all_gamma)
  #p <- ncol(all_gamma)
  
  # p x p covariance matrix
  Sigma <- covvar
  
  # Eigen decomposition (works for singular matrices)
  ee <- eigen(Sigma, symmetric = TRUE)
  lam <- ee$values
  vec <- ee$vectors
  
  # Keep only positive eigenvalues (numerically)
  keep <- lam > tol * max(lam)
  if (!any(keep)) stop("Covariance is numerically zero; all eigenvalues ~ 0.")
  
  lam_pos <- lam[keep]
  vec_pos <- vec[, keep, drop = FALSE]  # p x r
  r <- length(lam_pos)
  
  # Low-rank factor F: p x r, such that Sigma = F %*% t(F)
  # F = V * sqrt(Lambda)
  F <- sweep(vec_pos, 2, sqrt(lam_pos), `*`)
  
  # Now simulate: G: (nrep x r) with iid N(0,1)
  G <- matrix(rnorm(nrep * r), nrow = nrep, ncol = r)
  
  # Each row of V has distribution N(0, Sigma), possibly singular
  V <- G %*% t(F)   # (nrep x p)
  
  # Infinity norm per row
  Rfast::rowMaxs(abs(V), value = TRUE)
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
