library(tidyverse)
library(lightgbm)
library(vars) # technically unnessecary to impoort, but must be installed.

### source helper functions
source("helper_functions.R")
source("estimate_functions.R")



#todo: make sure packages are installed
# add docstrings and explanations
# consider whether or not remainders and such things should only be called if CCFs are true

estimate_stat <- function(data, L, B, 
                          # simulate mu and nu 
                          mu = matrix(rnorm(B), ncol = 1), #probably just as input
                          nu = matrix(rnorm(B * d), ncol = d), #prob just as input
                          num_rounds_for_train = 300,
                          # p = floor(log(n * L)),
                          lgb_params = list(),
                          objective = "regression",
                          parametric_plugin_AR1 = FALSE,
                          A = NULL,
                          remainder_true_ccfs = list(true_phi = NULL, true_psi = NULL)
                          # for AR(1) processes: 
                          # true ccf for phi should take only two arguments: a vector x, and a vector u, and a matrix A
                          # true ccf for psi should also take time input
                          # and return a length(x) times length(u) matrix
                          ) {
  #convert input data to a data.frame
  #allows for multiple input types, in particular data of type ts - time-series
  
  # browser()
  if(is.data.frame(data)) data <- data 
  else data <- data.frame(data)
  
  X <- data$X
  Y <- get_Y_mat(data)
  Ymat <- if (is.matrix(Y)) Y else cbind(Y)
  
  Tlen <- nrow(data)
  n <- Tlen / L
  d <- ncol(Ymat)
  
  p <- floor(log(n * L))
  
  Gamma_hat <- true_Gamma <- Gamma_parametric_plugin <- R1 <- R2 <- R3 <- complex(length.out = B)
  
  Covvar_Est <- list()
  all_gamma_hat <- matrix(NA, nrow = Tlen-L, ncol = 2*B)
  Covvar_Est_true <- list()
  all_gamma_true <- matrix(NA, nrow = Tlen-L, ncol = 2*B)
  Covvar_Est_parametric_plugin <- list()
  all_gamma_parametric_plugin <- matrix(NA, nrow = Tlen-L, ncol = 2*B)
  
  for (l in 1:L) {
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
    
    all_gamma_hat[((l - 1)*(n-1) +1) : (l*(n-1)), ] <- cbind(Re(lambda_hat), Im(lambda_hat))
    
    
    #estimate variance contribution for each l
    Covvar_Est[[l]] <- crossprod(cbind(Re(lambda_hat), Im(lambda_hat)))  # 2B x 2B
    
    Gamma_hat <- Gamma_hat + colSums(resX * resY)
    
    ### if true CCFs are give
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      #browser()
      
      
      true_phi <- remainder_true_ccfs$true_phi(X_eval_phi, mu, A)
      true_psi <- remainder_true_ccfs$true_psi(X_eval_psi, nu, A, index_eval_psi)
      
      true_resX <- cc_X - true_phi
      true_resY <- cc_Y - true_psi
      
      lambda_true <- true_resX * true_resY
      
      R1 <- R1 + colSums((cc_X - true_phi) * (true_psi - psi_hat_mat)) #psi residual
      R2 <- R2 + colSums((cc_Y - true_psi) * (true_phi - phi_hat_mat)) #phi residual
      R3 <- R3 + colSums((true_phi - phi_hat_mat) * (true_psi - psi_hat_mat)) #DML led
      
      true_Gamma <- true_Gamma + colSums((cc_X - true_phi) * (cc_Y - true_psi))
      
      all_gamma_true[((l - 1)*(n-1) +1) : (l*(n-1)), ] <- cbind(Re(lambda_true), Im(lambda_true))
      
      Covvar_Est_true[[l]] <- crossprod(cbind(Re(lambda_true), Im(lambda_true)))  # 2B x 2B
      
    }
    
    if(parametric_plugin_AR1) {
      
      #browser()
      
      fit <- vars::VAR(data)
      A_parametric_plugin <- Acoef(fit)[[1]]
      
      parametric_plugin_phi <- function(x, mu) remainder_true_ccfs$true_phi(x, mu, A_parametric_plugin)
      parametric_plugin_psi <- function(x, mu, t) remainder_true_ccfs$true_psi(x, mu, A_parametric_plugin, t)
      
      
      phi_parametric_plugin <- parametric_plugin_phi(X_eval_phi, mu)
      psi_parametric_plugin <- parametric_plugin_psi(X_eval_psi, nu, index_eval_psi)
      
      
      resX_parametric_plugin <- cc_X - phi_parametric_plugin
      resY_parametric_plugin <- cc_Y - psi_parametric_plugin
      
      
      lambda_parametric_plugin <- resX_parametric_plugin * resY_parametric_plugin
      
      Gamma_parametric_plugin <- Gamma_parametric_plugin + colSums(lambda_parametric_plugin)
      
      all_gamma_parametric_plugin[((l - 1)*(n-1) +1) : (l*(n-1)), ] <- cbind(Re(lambda_parametric_plugin), Im(lambda_parametric_plugin))
      
      Covvar_Est_parametric_plugin[[l]] <- crossprod(cbind(Re(lambda_parametric_plugin), Im(lambda_parametric_plugin)))  # 2B x 2B
      # could also return S_parametric_plugin if desired      
      
      
    }
    
    
  }
  
  normalizer <- (n - 1) * L
  
  Covvar_Est <- Reduce("+", Covvar_Est) / normalizer
  
  Gamma_hat <- Gamma_hat / normalizer
  R1 <- R1 / normalizer
  R2 <- R2 / normalizer
  R3 <- R3 / normalizer
  
  S_hat <- sqrt(normalizer) * max(abs(c(Re(Gamma_hat), Im(Gamma_hat))))
  
  
  output <- list(S_hat = S_hat, 
                 Covvar_Est = Covvar_Est,
                 all_gamma_hat = all_gamma_hat,
                 Gamma_hat = Gamma_hat
  )
  
  if(!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
    
    Covvar_Est_true <- Reduce("+", Covvar_Est_true) / normalizer
    
    true_Gamma <- true_Gamma / normalizer
    
    S_true <- sqrt(normalizer) * max(abs(c(Re(true_Gamma), Im(true_Gamma))))
    
    #normalize remainders similar to S
    R1 <- sqrt(normalizer) * max(abs(c(Re(R1), Im(R1))))
    R2 <- sqrt(normalizer) * max(abs(c(Re(R2), Im(R2))))
    R3 <- sqrt(normalizer) * max(abs(c(Re(R3), Im(R3))))
    
    output <- append(output,
                     list(Remainders = data.frame(R1 = R1, R2 = R2, R3 = R3), #used to be list
                          Covvar_Est_true = Covvar_Est_true,
                          all_gamma_true = all_gamma_true,
                          true_Gamma = true_Gamma,
                          S_true = S_true))
    }
  
  if(parametric_plugin_AR1) {
    Covvar_Est_parametric_plugin <- Reduce("+", Covvar_Est_parametric_plugin) / normalizer
    
    Gamma_parametric_plugin <- Gamma_parametric_plugin / normalizer
    
    S_parametric_plugin <- sqrt(normalizer) * max(abs(c(Re(Gamma_parametric_plugin), Im(Gamma_parametric_plugin))))
    
    output <- append(output,
                     list(parametric_plugin = 
                            list(Covvar_Est_parametric_plugin = Covvar_Est_parametric_plugin,
                                 Gamma_parametric_plugin = Gamma_parametric_plugin,
                                 all_gamma_parametric_plugin = all_gamma_parametric_plugin,
                                 S_parametric_plugin = S_parametric_plugin)
                          )
                     )
    
  }
  
  return(output)
}


# data_test <- simulate_AR_process(100, A)
# test <- estimate_stat(data = data_test, L = 10, B = 2
#               # remainder_true_ccfs = list(
#               #   true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#               #   true_psi = function(x, u, A, t) char_func_cond_Y_given_X_mat(A, t, x, u)
#               # )
#             )
# 
# 
# small_data <- simulate_AR_process(10, d = 3, A = matrix(c(0.1, rep(0, 3), rep(0.1, 12)), nrow = 4, byrow = T))
# estimate_stat(small_data, 5, 2)
# 
# estimate_stat(my_X$discretized_path, 2, 2)





