
library(tidyverse)
library(lightgbm)
library(vars) # technically unnessecary to impoort, but must be installed.

### source helper functions
source("helper_functions.R")
source("estimate_functions.R")



#todo: make sure packages are installed
# add docstrings and explanations
# consider whether or not remainders and such things should only be called if CCFs are true

estimate_stat_oracle_parametric_only <- function(data, L, B, 
                          # simulate mu and nu 
                          mu = matrix(rnorm(B), ncol = 1), #probably just as input
                          nu = matrix(rnorm(B * d), ncol = ncol(data)-1), #prob just as input
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
  if(!is.data.frame(data)) data <- as.data.frame(data)
  
  X <- data$X
  Y <- get_Y_mat(data)
  Ymat <- if (is.matrix(Y)) Y else cbind(Y)
  
  Tlen <- nrow(data)
  n <- Tlen / L
  d <- ncol(Ymat)
  
  p <- log(n * L)
  
  true_Gamma <- Gamma_parametric_plugin <- complex(length.out = B)
  
  Covvar_Est_true <- list()
  Covvar_Est_parametric_plugin <- list()

  
  fit <- vars::VAR(data)
  A_parametric_plugin <- Acoef(fit)[[1]]
  
  stationary_covariance <- stationary_covariance(A_parametric_plugin)
  
  true_phi_func <- function(x, mu) remainder_true_ccfs$true_phi(x, mu, A)
  true_psi_func <- function(x, mu, t) remainder_true_ccfs$true_psi(x, mu, A, t, stationary_covariance = stationary_covariance)
  
  parametric_plugin_phi <- function(x, mu) remainder_true_ccfs$true_phi(x, mu, A_parametric_plugin)
  parametric_plugin_psi <- function(x, mu, t) remainder_true_ccfs$true_psi(x, mu, A_parametric_plugin, t, stationary_covariance = stationary_covariance)
  
  
  
  for (l in 1:(L - 1)) {
    
    
    #browser()
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
    
    ### empirical CF on evaluation for X
    cc_X <- matrix(
      exp(1i * rep(mu, each = N_eval_phi) * rep(X_shifted_eval_phi, times = B)),
      nrow = N_eval_phi, ncol = B)
    
    ### stuff for psi
    #X_train_psi <- X[index_train]
    #Y_train_psi <- Ymat[index_train, , drop = F] 
    
    index_eval_psi <- index_eval[1:(n - 1)]
    
    X_eval_psi <- X[index_eval_psi]
    Y_eval_psi <- Ymat[index_eval_psi, , drop = F] 
    
    ### empirical CF on evaluation for Y
    cc_Y <- exp(1i * (Y_eval_psi %*% t(nu)))
 
    ### if true CCFs are give
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      #browser()
      
      
      true_phi <- true_phi_func(X_eval_phi, mu)
      true_psi <- true_psi_func(X_eval_psi, nu, index_eval_psi)
      
      true_resX <- cc_X - true_phi
      true_resY <- cc_Y - true_psi
      
      lambda_true <- true_resX * true_resY
      
      true_Gamma <- true_Gamma + colSums((cc_X - true_phi) * (cc_Y - true_psi))
      
      Covvar_Est_true[[l]] <- crossprod(cbind(Re(lambda_true), Im(lambda_true)))  # 2B x 2B
      
    }
    
    if(parametric_plugin_AR1) {
    
      
      phi_parametric_plugin <- parametric_plugin_phi(X_eval_phi, mu)
      psi_parametric_plugin <- parametric_plugin_psi(X_eval_psi, nu, index_eval_psi)
      
      
      resX_parametric_plugin <- cc_X - phi_parametric_plugin
      resY_parametric_plugin <- cc_Y - psi_parametric_plugin
      
      
      lambda_parametric_plugin <- resX_parametric_plugin * resY_parametric_plugin
      
      Gamma_parametric_plugin <- Gamma_parametric_plugin + colSums(lambda_parametric_plugin)
      
      
      Covvar_Est_parametric_plugin[[l]] <- crossprod(cbind(Re(lambda_parametric_plugin), Im(lambda_parametric_plugin)))  # 2B x 2B
    
    }
    
    
  }
  
  normalizer <- (n - 1) * (L - 1)
  
  # Covvar_Est <- Reduce("+", Covvar_Est) / normalizer
  # 
  # Gamma_hat <- Gamma_hat / normalizer
  # R1 <- R1 / normalizer
  # R2 <- R2 / normalizer
  # R3 <- R3 / normalizer
  # 
  # S_hat <- sqrt(normalizer) * max(abs(c(Re(Gamma_hat), Im(Gamma_hat))))
  # 
  
  output <- list(#S_hat = S_hat, 
                 #Covvar_Est = Covvar_Est,
                 #Gamma_hat = Gamma_hat
                 )
  
  if(!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
    
    Covvar_Est_true <- Reduce("+", Covvar_Est_true) / normalizer
    
    true_Gamma <- true_Gamma / normalizer
    
    S_true <- sqrt(normalizer) * max(abs(c(Re(true_Gamma), Im(true_Gamma))))
    
    # #normalize remainders similar to S
    # R1 <- sqrt(normalizer) * max(abs(c(Re(R1), Im(R1))))
    # R2 <- sqrt(normalizer) * max(abs(c(Re(R2), Im(R2))))
    # R3 <- sqrt(normalizer) * max(abs(c(Re(R3), Im(R3))))
    
    output <- append(output,
                     list(#Remainders = data.frame(R1 = R1, R2 = R2, R3 = R3), #used to be list
                          Covvar_Est_true = Covvar_Est_true,
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
                                 S_parametric_plugin = S_parametric_plugin)
                          )
                     )
    
  }
  
  #browser()
  
  return(output)
}
# 
# source("simulate_AR_process.R")
# source("conditional_distributions.R")
# 
# 
# set.seed(420)
# A_mat <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)
# 
# 
# data <- simulate_AR_process(1e3, d = 3, A = A_mat)
# 
# 
# 
# estimate_stat_oracle_parametric_only(data, n, L = 10, B = 2,
#               A = A_mat,
#                mu = matrix(rnorm(2), ncol = 1), #probably just as input
#                nu = matrix(rnorm(2 * 3), ncol = 3),
#               remainder_true_ccfs = list(
#                 true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                 true_psi = function(x, u, A, t) char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#               )
#             )

# source("estimate_test.R")
# 
# 
# estimate_stat_oracle_parametric_only(
#   data = data, L = 10, B = 2,
#   mu =  matrix(rnorm(2), ncol = 1), #probably just as input
#   nu = matrix(rnorm(2 * 3), ncol = 3),
#   A = A_mat,
#   parametric_plugin_AR1 = TRUE,
#   remainder_true_ccfs = list(
#     true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#     true_psi = function(x, u, A, t) char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#   )
# )
# 
# 

# estimate_stat(small_data, 5, 2)

# estimate_stat(my_X$discretized_path, 2, 2)





