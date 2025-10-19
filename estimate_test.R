library(tidyverse)
library(lightgbm)

### source helper functions
source("helper_functions.R")
source("estimate_functions.R")

estimate_stat <- function(data, n, L, B, 
                          # simulate mu and nu 
                          mu = matrix(rnorm(B), ncol = 1), #probably just as input
                          nu = matrix(rnorm(B * d), ncol = d), #prob just as input
                          num_rounds_for_train = 300,
                          p = 5,
                          lgb_params = list(),
                          objective = "regression") {
  # browser()
  X <- data$X
  Y <- data$Y #perhaps use the function that greps Y-matrix
  Ymat <- if (is.matrix(Y)) Y else cbind(Y)
  
  Tlen <- nrow(data)
  d <- ncol(Ymat)
  
  Gamma <- complex(length.out = B)
  Covvar_Est <- list()
  
  for (l in 1:(L - 1)) {
    # indices for training and evaluation
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
    
    phi_hat <- phi_hat(X_train_phi, 
                       X_shifted_train_phi,
                       X_eval_phi, 
                       mu, 
                       num_rounds_for_train,
                       lgb_params,
                       objective)
    
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
    
    psi_hat <- psi_hat(X_train_psi, 
                       Y_train_psi, 
                       X_eval_psi, 
                       nu, 
                       num_rounds_for_train, 
                       lgb_params,
                       objective)
    
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
  
  normalizer <- (n - 1) * (L - 1)
  
  Covvar_Est <- Reduce("+", Covvar_Est) / normalizer
  
  Gamma_hat <- Gamma / normalizer
  S_hat <- sqrt(normalizer) * max(abs(c(Re(Gamma_hat), Im(Gamma_hat))))
  
  return(list( S_hat = S_hat, Covvar_Est = Covvar_Est))
}






