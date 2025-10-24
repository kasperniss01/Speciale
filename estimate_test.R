library(tidyverse)
library(lightgbm)

### source helper functions
source("helper_functions.R")
source("estimate_functions.R")

estimate_stat <- function(data, L, B, 
                          # simulate mu and nu 
                          mu = matrix(rnorm(B), ncol = 1), #probably just as input
                          nu = matrix(rnorm(B * d), ncol = d), #prob just as input
                          num_rounds_for_train = 300,
                          # p = floor(log(n * L)),
                          lgb_params = list(),
                          objective = "regression",
                          remainder_true_ccfs = list(true_phi = NULL, true_psi = NULL)
                          # for AR(1) processes: 
                          # true ccf for phi should take only two arguments: a vector x, and a vector u
                          # true ccf for psi should also take time input
                          # and return a length(x) times length(u) matrix
                          ) {
  # browser()
  X <- data$X
  Y <- get_Y_mat(data)
  Ymat <- if (is.matrix(Y)) Y else cbind(Y)
  
  Tlen <- nrow(data)
  n <- Tlen / L
  d <- ncol(Ymat)
  
  p <- log(n * L)
  
  Gamma_hat <- R1 <- R2 <- R3 <- true_Gamma <- complex(length.out = B)
  
  Covvar_Est <- list()
  Covvar_Est_true <- list()
  
  for (l in 1:(L - 1)) {
    # indices for training and evaluation
    index_train <- II_bar_l(Tlen, L, l, p)
    index_eval <- II_l(Tlen, L, l)
    
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
    
    phi_hat_mat <- phi_hat(X_train_phi, 
                       X_shifted_train_phi,
                       X_eval_phi, 
                       mu, 
                       num_rounds_for_train,
                       lgb_params,
                       objective)
    
    ### empirical CF on evaluation for X
    cc_X <- matrix(
      exp(1i * rep(mu, each = N_eval_phi) * rep(X_shifted_eval_phi, times = B)),
      nrow = N_eval_phi, ncol = B)
    
    ### stuff for psi
    X_train_psi <- X[index_train]
    Y_train_psi <- Ymat[index_train, , drop = F] 
    
    index_eval_psi <- index_eval[1:(n - 1)]
    
    X_eval_psi <- X[index_eval_psi]
    Y_eval_psi <- Ymat[index_eval_psi, , drop = F] 
    
    psi_hat_mat <- psi_hat(X_train_psi, 
                       Y_train_psi, 
                       X_eval_psi, 
                       nu, 
                       num_rounds_for_train, 
                       lgb_params,
                       objective)
    
    ### empirical CF on evaluation for Y
    cc_Y <- exp(1i * (Y_eval_psi %*% t(nu)))
    
    ### calculate residual terms
    resX <- cc_X - phi_hat_mat
    resY <- cc_Y - psi_hat_mat
    
    # browser()
    lambda_hat <- resX * resY
    
    #estimate variance contribution for each l
    Covvar_Est[[l]] <- crossprod(cbind(Re(lambda_hat), Im(lambda_hat)))  # 2B x 2B
    
    Gamma_hat <- Gamma_hat + colSums(resX * resY)
    
    ### if true CCFs are give
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      # browser()
      
      true_phi <- remainder_true_ccfs$true_phi(X_eval_phi, mu)
      true_psi <- remainder_true_ccfs$true_psi(X_eval_psi, nu, index_eval_psi)
      
      true_resX <- cc_X - true_phi
      true_resY <- cc_Y - true_psi
      
      lambda_true <- true_resX * true_resY
      
      R1 <- R1 + colSums((true_phi - phi_hat_mat) * (true_psi - psi_hat_mat))
      R2 <- R2 + colSums((cc_X - true_phi) * (true_psi - psi_hat_mat))
      R3 <- R3 + colSums((cc_Y - true_psi) * (true_phi - phi_hat_mat))
      
      true_Gamma <- true_Gamma + colSums((cc_X - true_phi) * (cc_Y - true_psi))
      
      Covvar_Est_true[[l]] <- crossprod(cbind(Re(lambda_true), Im(lambda_true)))  # 2B x 2B
      
    }
    
  }
  
  normalizer <- (n - 1) * (L - 1)
  
  Covvar_Est <- Reduce("+", Covvar_Est) / normalizer
  
  Gamma_hat <- Gamma_hat / normalizer
  R1 <- R1 / normalizer
  R2 <- R2 / normalizer
  R3 <- R3 / normalizer
  
  S_hat <- sqrt(normalizer) * max(abs(c(Re(Gamma_hat), Im(Gamma_hat))))
  
  
  output <- list(S_hat = S_hat, 
                 Covvar_Est = Covvar_Est,
                 Gamma_hat = Gamma_hat
                 )
  
  if(!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
    
    Covvar_Est_true <- Reduce("+", Covvar_Est_true) / normalizer
    
    true_Gamma <- true_Gamma / normalizer
    
    S_true <- sqrt(normalizer) * max(abs(c(Re(true_Gamma), Im(true_Gamma))))
    
    output <- append(output,
                     list(Remainders = list(R1 = R1, R2 = R2, R3 = R3),
                          Covvar_Est_true = Covvar_Est_true,
                          true_Gamma = true_Gamma,
                          S_true = S_true))
  }
  
  return(output)
}

# test <- estimate_stat(data2, n, L, B,
#               remainder_true_ccfs = list(
#                 true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A[1,1], x, u),
#                 true_psi = function(x, u, t) char_func_cond_Y_given_X_mat(A[1,1], A[2,1], A[2,2], t, x, u)
#               )
#             )


# small_data <- simulate_AR_process(10, d = 3, A = matrix(c(0.1, rep(0, 3), rep(0.1, 12)), nrow = 4, byrow = T))
# estimate_stat(small_data, 5, 2)





