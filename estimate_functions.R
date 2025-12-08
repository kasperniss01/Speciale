### ------ functions for estimation of phi and psi -----

### todo: add docstrings 

### fitting a lightgbm regression model that defaults to squared-error-loss ###
fit_lgb <- function(X, y, num_round = 150, params = list(),
                    objective = "regression") {
  
  # X covariates
  # y response
  # objective is what kind of regression to do
  
  X_mat <- as.matrix(X)
  n <- nrow(X_mat) # number of observations
  
  dtrain <- lgb.Dataset(data = as.matrix(X), label = y)
  
  default <- list(
    objective = objective,
    learning_rate = 0.05,  #hopefully sane parameters 
    num_leaves = 20,
    min_data_in_leaf = max(100, floor(0.005 * n)),
    feature_fraction = 1,
    bagging_fraction = 0.7,
    bagging_freq = 1,
    lambda_l2 = 2,
    verbose = -1,
    seed = 1
  )
  
  params <- modifyList(default, params)
  lgb.train(params = params, data = dtrain, nrounds = num_round) #fits lgb model
}


### ------------ estimates of phi and psi -------------- ###

### phi_hat trained on training_X, evaluated on evaluation_X for multiple mu ###

phi_hat <- function(training_X, training_X_shifted , evaluation_X, mu, num_rounds, lgb_params,
                    objective = "regression") {
  # training X are X to be trained on
  # evaluation X are X to be evaluated on
  # mu are evaluation points for characteristic functions 
    # also used as covariates!
  
  # returns matrix of predictions based on lgbm-model
  
  N_train <- length(training_X)
  N_eval <- length(evaluation_X)
  B <- if (is.matrix(mu)) nrow(mu) else length(mu)
  
  mu_train <- mu[rep(seq_len(B), each = N_train), , drop = FALSE]
  colnames(mu_train) <- paste0("mu", seq_len(1))
  
  features_train <- cbind(
    x = rep(training_X, times = B),
    mu = mu_train
  )
  
  mu_eval <- mu[rep(seq_len(B), each = N_eval), , drop = FALSE]
  colnames(mu_eval) <- paste0("mu", seq_len(1))
  
  features_eval <- cbind(
    x = rep(evaluation_X, times = B),
    mu = mu_eval
  )
  
  real_response <- cos(rep(mu, each = N_train) * rep(training_X_shifted, times = B))
  imag_response <- sin(rep(mu, each = N_train) * rep(training_X_shifted, times = B))
  
  #fit models for real and imaginary response
  model_real <- fit_lgb(features_train, 
                        real_response,
                        num_rounds,
                        lgb_params,
                        objective)
  model_imag <- fit_lgb(features_train,
                        imag_response,
                        num_rounds,
                        lgb_params,
                        objective)
  
  real_hat <- predict(model_real, features_eval)
  imag_hat <- predict(model_imag, features_eval)
  
  #output predictions
  matrix(real_hat, nrow = N_eval, ncol = B) + 1i * matrix(imag_hat, nrow = N_eval, ncol = B)
}

### function to create estimate of psi across multiple values of nu
psi_hat <- function(training_X, training_Y , evaluation_X, nu, num_rounds, lgb_params,
                    objective = "regression") {
  # training X and Y are X and Y to be trained on
  # evaluation X are X to be evaluated on
  # nu are evaluation points for characteristic functions 
  # also used as covariates!
  
  # returns matrix of predictions based on lgbm-model
  
  # browser()
  
  N_train <- length(training_X)
  N_eval <- length(evaluation_X)
  B <- if (is.matrix(nu)) nrow(nu) else length(nu)
  d <- ncol(training_Y)
  
  nu_train <- nu[rep(seq_len(B), each = N_train), , drop = FALSE]
  colnames(nu_train) <- paste0("nu", seq_len(d))
  
  features_train <- cbind(
    x = rep(training_X, times = B),
    nu_train
  )
  
  nu_mult_Y_train <-training_Y %*% t(nu)
  
  real_response <- cos(as.vector(nu_mult_Y_train))
  imag_response <- sin(as.vector(nu_mult_Y_train))
  
  #fit models for real and imaginary response
  model_real <- fit_lgb(features_train,
                        real_response,
                        num_rounds,
                        params = lgb_params,
                        objective)
  model_imag <- fit_lgb(features_train,
                        imag_response,
                        num_rounds,
                        params = lgb_params,
                        objective)
  
  # rep nu for evaluation
  nu_eval <- nu[rep(seq_len(B), each = N_eval), , drop = FALSE]
  colnames(nu_eval) <- paste0("nu", seq_len(d))
  
  features_eval <- cbind(
    x = rep(evaluation_X, times = B),
    nu_eval
  )
  
  # perform predictions
  real_hat <- predict(model_real, as.matrix(features_eval))
  imag_hat <- predict(model_imag, as.matrix(features_eval))
  
  #ouput predictions
  matrix(real_hat, nrow = N_eval, ncol = B) + 1i * matrix(imag_hat, nrow = N_eval, ncol = B)
  
}
