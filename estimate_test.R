library(tidyverse)
library(lightgbm)
rm(list = ls(
))
setwd("/Users/davidsenderovitz/Library/CloudStorage/OneDrive-Personligt/Dokumenter/Universitet/Speciale/Speciale_git/Speciale")

source("simulate_AR_process.R")

### function to split timeseries into L blocks of size n
make_blocks <- function(Tlen, n, L) {
  stopifnot(Tlen == n * L)
  split(seq_len(Tlen), rep(seq_len(L), each = n))
}

### fetch ell'th (l) block from the L blocks - corresponds to evaluation data
II_l <- function(Tlen, n, L, l) { 
  stopifnot(1 <= l & l <= L)
  make_blocks(Tlen, n, L)[[l]]
}

### fetch l'th training data according to equation ??
II_bar_l <- function(Tlen, n, L, l, p) {
  idx_eval <- II_l(Tlen, n, L, l)
  
  left  <- max(1, min(idx_eval) - p)
  right <- min(Tlen, max(idx_eval) + p)
  
  setdiff(seq_len(Tlen), left:right)
}

#get pairs (t, t + 1) such that both entries are training points
pairs_from_mask <- function(in_mask) {
  Tlen <- length(in_mask)
  which(in_mask[-Tlen] & in_mask[-1])
}

# lightgbm wrapper with sane defaults
fit_lgb <- function(X, y, num_round = 300, params = list()) {
  dtrain <- lgb.Dataset(data = as.matrix(X), label = y)
  default <- list(
    objective = "regression_l2",
    learning_rate = 0.05, # what is all this
    num_leaves = 31,
    min_data_in_leaf = 20,
    feature_fraction = 0.9,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    verbose = -1,
    seed = 1
  )
  params <- modifyList(default, params)
  lgb.train(params = params, data = dtrain, nrounds = num_round)
}

estimate_stat <- function(data, n, L, B, 
                          num_rounds_for_train = 300,
                          lgb_params = list()) {
  # browser()
  X <- data$X
  Y <- data$Y #perhaps use the function that greps Y-matrix
  Ymat <- if (is.matrix(Y)) Y else cbind(Y)
  
  Tlen <- nrow(data)
  d <- ncol(Ymat)
  
  # simulate mu and nu 
  mu <- rnorm(B)
  nu <- matrix(rnorm(B * d), ncol = d)
  
  I <- complex(real = 0, imaginary = 1)
  
  Gamma <- complex(length.out = B)
  Covvar_est <- matrix(0, 2 * B, 2 * B)
  
  
  p <- 5
  
  for (l in 1:(L - 1)) {
    browser()
    
    # indices for training and evalution
    index_train <- II_bar_l(Tlen, n, L, l, p)
    index_eval <- II_l(Tlen, n, L, l)
    
    #masking indices
    in_train <- rep(FALSE, Tlen)
    in_train[index_train] <- TRUE
    
    ### φ training
    # φ training pairs: both t and t+1 must be in training region
    index_train_phi  <- pairs_from_mask(in_train) 
    shifted_index_train_phi <- index_train_phi + 1
    
    X_train_phi <- X[index_train_phi]
    X_shifted_train_phi <- X[shifted_index_train_phi]
    
    N_train_phi <- length(X_train_phi)
    
    phi_features_train <- cbind(
      x = rep(X_train_phi, times = B),
      mu = rep(mu, each = N_train_phi)
    )
    
    phi_real_response <- cos(rep(mu, each = N_train_phi) * rep(X_shifted_train_phi, times = B))
    phi_imag_response <- sin(rep(mu, each = N_train_phi) * rep(X_shifted_train_phi, times = B))
    
    model_phi_real <- fit_lgb(phi_features_train, phi_real_response, num_round = num_rounds_for_train, params = lgb_params)
    model_phi_imag <- fit_lgb(phi_features_train, phi_imag_response, num_round = num_rounds_for_train, params = lgb_params)
    
    ### φ predictions
    index_eval_phi <- index_eval[1:(n - 1)]
    shifted_index_eval_phi <- index_eval_phi + 1
    
    X_eval_phi <- X[index_eval_phi]
    X_shifted_eval_phi <- X[shifted_index_eval_phi]
    
    N_eval_phi <- length(index_eval_phi)
    
    phi_features_eval <- cbind(
      x = rep(X_eval_phi, times = B),
      mu = rep(mu, each = N_eval_phi)
    )
    
    # perform predictions
    phi_real_hat_vec <- predict(model_phi_real, as.matrix(phi_features_eval))
    phi_imag_hat_vec <- predict(model_phi_imag, as.matrix(phi_features_eval))
    
    #combine into complex estimates
    phi_hat <- matrix(phi_real_hat_vec, nrow = N_eval_phi, ncol = B) +
      I * matrix(phi_imag_hat_vec, nrow = N_eval_phi, ncol = B)
    
    ### true CCF on evaluation for X
    cc_X <- matrix(
      exp(I * rep(mu, each = N_eval_phi) * rep(X_shifted_eval_phi, times = B)),
      nrow = N_eval_phi, ncol = B
    )
    
    ### ψ training
    # ψ training pairs: both X and Y in training
    X_train_psi <- X[index_train]
    Y_train_psi <- Ymat[index_train, , drop = F] 
    
    N_train_psi <- length(X_train_psi)
    
    # Expand ν into d columns with consistent names
    nu_train_rep <- nu[rep(seq_len(B), each = N_train_psi), , drop = FALSE]
    colnames(nu_train_rep) <- paste0("nu", seq_len(d))
    
    psi_features_train <- cbind(
      x = rep(X_train_psi, times = B),
      nu_train_rep
    )
    
    nu_mult_Y_train <-Y_train_psi %*% t(nu)
    
    psi_real_response <- cos(as.vector(nu_mult_Y_train))
    psi_imag_response <- sin(as.vector(nu_mult_Y_train))
    
    model_psi_real <- fit_lgb(psi_features_train, 
                              psi_real_response, 
                              num_round = num_rounds_for_train, 
                              params = lgb_params)
    model_psi_imag <- fit_lgb(psi_features_train, 
                              psi_imag_response, 
                              num_round = num_rounds_for_train, 
                              params = lgb_params)
    
    ### ψ predictions
    index_eval_psi <- index_eval[1:(n - 1)]
    shifted_index_eval_psi <- index_eval_psi + 1
    
    X_eval_psi <- X[index_eval_psi]
    Y_eval_psi <- Ymat[index_eval_psi, , drop = F] 
    
    N_eval_psi <- length(index_eval_psi)
    
    # rep nu for evaluation
    nu_eval_rep <- nu[rep(seq_len(B), each = N_eval_psi), , drop = FALSE]
    colnames(nu_eval_rep) <- paste0("nu", seq_len(d))
    
    psi_features_eval <- cbind(
      x = rep(X_eval_psi, times = B),
      nu_eval_rep
    )
    
    # perform predictions
    psi_real_hat_vec <- predict(model_psi_real, as.matrix(psi_features_eval))
    psi_imag_hat_vec <- predict(model_psi_imag, as.matrix(psi_features_eval))
    
    #combine into complex estimates
    psi_hat <- matrix(psi_real_hat_vec, nrow = N_eval_psi, ncol = B) +
      I * matrix(psi_imag_hat_vec, nrow = N_eval_psi, ncol = B)
    
    ### true CCF on evaluation for Y
    cc_Y <- exp(I * (Y_eval_psi %*% t(nu)))
    
    ### calculate residual terms
    resX <- cc_X - phi_hat
    resY <- cc_Y - psi_hat
    
    lambda <- resX * resY
    
    Lambda_mat <- cbind(Re(lambda), Im(lambda))
    var_contrib <- crossprod(Lambda_mat)
    
    Covvar_est <- Covvar_est +  var_contrib
    
    
    Gamma <- Gamma + colSums(resX * resY)
  }
  
  Covvar_est <- Covvar_est / ((n - 1) * (L - 1))
  Gamma_hat <- Gamma / ((n - 1) * (L - 1))
  
  
  S_hat <- sqrt((n - 1) * (L - 1)) * max(abs(Re(Gamma_hat)), abs(Im(Gamma_hat)))
  
  return(list( S_hat = S_hat, Covvar_est = Covvar_est))
}



est_stat <- estimate_stat(data = sim, n = 10, L = 3, B = 5, lgb_params = list(learning_rate = 0.05,
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
n_sim <- 10
S_hats <- numeric(n_sim)


reject_values <- numeric(n_sim)
for(k in 1:n_sim) {
  data_sim <- simulate_AR_process(gamma = 0, 
                                  n = Tlen, 
                                  burnin = 0, 
                                  A = matrix(c(0.1, 0.7, -0.1, 0.3), nrow = 2))
  
  est_stat_sim <- estimate_stat(data = data_sim, n = 100, L = 5, B = 10)
  S_hat_sim <- est_stat_sim$S_hat
  Covvar_est_sim <- est_stat_sim$Covvar_est
  
  crit_val <- sim_crit_value(covvar_est = Covvar_est_sim, alpha = 0.05)
  
  reject_values[k] <- S_hat_sim > crit_val
  S_hats[k] <- S_hat_sim
  if(k %% 10 == 0) print(k)


}

reject_values %>% mean

ggplot() +
  geom_histogram(aes(S_hats), bins = 30, color = "white")





ggplot() +
  geom_histogram(aes(
    rnorm(20*1000, sd = 0.1) %>% 
      matrix(ncol = 20) %>% 
      abs() %>% 
      apply(1, max)), bins = 30, color = "white")









