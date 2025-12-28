### test the hypothesis ###

source("estimate_test.R")
source("sim_crit_value.R")
source("helper_functions.R")

test_hypothesis <- function(data, L, B, mu, nu, alpha) {
  #if H_0 is to be rejected or not
  #input: data, L integer: # of folds, B is # of evaluation points
  # mu and nu are evaluation points
  #alpha is the desired significance level
  
  #output: list of class "hyp_test"
    # call, reject, p_val, S_hat, critical value and significance level
  est <- estimate_stat(
    data = data, 
    L = L, 
    B = B,
    mu = mu, 
    nu = nu
  )
  
  S_hat <- est$S_hat
  covvar_est <- est$Covvar_Est
  
  #simulate critical value
  draws_alt <- sim_crit_draws_alt2(covvar_est, nrep = 1e5)
  crits_alt <- crit_from_draws(draws_alt, alpha)
  
  p_val <- (sum(draws_alt >= S_hat) + 1) / (length(draws_alt) + 1)
  
  #reject or not 
  reject <- S_hat > crits_alt #p_val <= alpha - equivalent
  
  #outout
  out <- list(
    call = match.call(),
    reject = reject,
    p_val = p_val, 
    S_hat = S_hat,
    crit_value = crits_alt,
    alpha = alpha
  )
  
  class(out) <- "hyp_test"
  out
  
}

