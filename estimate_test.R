library(tidyverse)
library(lightgbm)

### source helper functions
source("helper_functions.R")

estimate_stat <- function(data, n, L, B, 
                          num_rounds_for_train = 300,
                          p = 5,
                          lgb_params = list()) {
  X <- data$X
  Y <- data$Y #perhaps use the function that greps Y-matrix
  Ymat <- if (is.matrix(Y)) Y else cbind(Y)
  
  Tlen <- nrow(data)
  d <- ncol(Ymat)
  
  # simulate mu and nu 
  mu <- matrix(rnorm(B), ncol = 1)
  nu <- matrix(rnorm(B * d), ncol = d)
  
  Gamma <- complex(length.out = B)
  Covvar_Est <- list()
  
  for (l in 1:(L - 1)) {
    # indices for training and evalution
    index_train <- II_bar_l(Tlen, n, L, l, p)
    index_eval <- II_l(Tlen, n, L, l)
    
    ### stuff for phi
    #masking indices
    in_train <- rep(FALSE, Tlen)
    in_train[index_train] <- TRUE
    
    # phi training pairs: both t and t+1 must be in training region
    index_train_phi  <- pairs_from_mask(in_train) 
    shifted_index_train_phi <- index_train_phi + 1
    
    X_train_phi <- X[index_train_phi]
    X_shifted_train_phi <- X[shifted_index_train_phi]
    
    index_eval_phi <- index_eval[1:(n - 1)]
    shifted_index_eval_phi <- index_eval_phi + 1
    
    X_eval_phi <- X[index_eval_phi]
    X_shifted_eval_phi <- X[shifted_index_eval_phi]
    
    N_eval_phi <- length(X_eval_phi)
    
    phi_hat <- phi_hat(X_train_phi, X_shifted_train_phi, X_eval_phi, mu, num_rounds_for_train, lgb_params)
    
    ### true CCF on evaluation for X
    cc_X <- matrix(
      exp(1i * rep(mu, each = N_eval_phi) * rep(X_shifted_eval_phi, times = B)),
      nrow = N_eval_phi, ncol = B)
    
    ### stuff for psi
    X_train_psi <- X[index_train]
    Y_train_psi <- Ymat[index_train, , drop = F] 
    
    index_eval_psi <- index_eval[1:(n - 1)]
    
    X_eval_psi <- X[index_eval_psi]
    Y_eval_psi <- Ymat[index_eval_psi, , drop = F] 
    
    psi_hat <- psi_hat(X_train_psi, Y_train_psi, X_eval_psi, nu, num_rounds_for_train, lgb_params)
    
    ### true CCF on evaluation for Y
    cc_Y <- exp(1i * (Y_eval_psi %*% t(nu)))
    
    ### calculate residual terms
    resX <- cc_X - phi_hat
    resY <- cc_Y - psi_hat
    
    lambda <- resX * resY
    
    var_contrib <- matrix(0, 2 * B, 2 * B)
    for(t in 1:(n - 1)) {
      vv <- c(Re(lambda[t, ]), Im(lambda[t, ]))
      var_contrib <- var_contrib + (vv %o% vv)

    }
    Covvar_Est[[l]] <- var_contrib
    
    Gamma <- Gamma + colSums(resX * resY)
  }
  
  Covvar_Est <- Reduce("+", Covvar_Est) / ((n - 1) * (L - 1))
  
  Gamma_hat <- Gamma / ((n - 1) * (L - 1))
  S_hat <- sqrt((n - 1)*(L - 1)) * max(abs(c(Re(Gamma_hat), Im(Gamma_hat))))
  
  return(list( S_hat = S_hat, Covvar_Est = Covvar_Est))
}



est_stat <- estimate_stat(data = sim, n = 100, L = 10, B = 5, lgb_params = list(learning_rate = 0.05,
                                                                  num_leaves = 31,
                                                                  verbose = -1))

sim_crit_value <- function(n = 1e3, covvar_est, alpha = 0.05) {
  B <- nrow(covvar_est) / 2
  dist <- replicate(n, {
    Z <- rnorm(2 * B)
    
    root_covvar_est <- t(chol(covvar_est))
    max(abs(root_covvar_est %*% Z))
  })
  
  unname(quantile(dist, 1 - alpha/2))
}

sim_crit_value(1e5, est_stat$Covvar_Est, 0.05)


sim <- simulate_AR_process(gamma = 0, 
                           n = 30, 
                           burnin = 50, 
                           A = matrix(c(0.1, 0.7, -0.1, 0.3), nrow = 2))

Tlen <- 500
rej_values <- replicate(2, {
  data_sim <- simulate_AR_process(gamma = 0, 
                                  n = Tlen, 
                                  burnin = 0, 
                                  A = matrix(c(0.1, 0.7, -0.1, 0.3), nrow = 2))
  
  est_stat_sim <- estimate_stat(data = data_sim, n = 100, L = 5, B = 10)
  S_hat_sim <- est_stat_sim$S_hat
  Covvar_est_sim <- est_stat_sim$Covvar_Est
  
  crit_val <- sim_crit_value(covvar_est = Covvar_est_sim, alpha = 0.05)
  
  reject <- S_hat_sim > crit_val

  
  list(reject = reject, S_hat = S_hat_sim)

}, simplify = F)


S_hat_vals <- c()
for (i in 1:100) {
  S_hat_vals[i] <- rej_values[[i]]$S_hat
}

ggplot() + 
  geom_histogram(aes(x = S_hat_vals, y = after_stat(density)), color = "white")




