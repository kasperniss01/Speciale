### -------- helper functions ------- 

### todo: update functions with docline and parameter input assumptions

## ------ functions for estimation of phi and psi -----

### fitting a lightgbm regression model med squared-error-loss  
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

### function to create estimate of phi across multiple values of mu
phi_hat <- function(training_X, training_X_shifted , evaluation_X, mu, num_rounds, lgb_params) {
  N_train <- length(training_X)
  N_eval <- length(evaluation_X)
  B <- nrow(mu)
  # browser()
  
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
  
  model_real <- fit_lgb(features_train, real_response, num_rounds, lgb_params)
  model_imag <- fit_lgb(features_train, imag_response, num_rounds, lgb_params)
  
  real_hat <- predict(model_real, features_eval)
  imag_hat <- predict(model_imag, features_eval)
  
  matrix(real_hat, nrow = N_eval, ncol = B) + 1i * matrix(imag_hat, nrow = N_eval, ncol = B)
}

### function to create estimate of psi across multiple values of nu
psi_hat <- function(training_X, training_Y , evaluation_X, nu, num_rounds, lgb_params) {
  N_train <- length(training_X)
  N_eval <- length(evaluation_X)
  B <- nrow(nu)
  d <- ncol(training_Y)
  # browser()
  
  nu_train <- nu[rep(seq_len(B), each = N_train), , drop = FALSE]
  colnames(nu_train) <- paste0("nu", seq_len(d))
  
  features_train <- cbind(
      x = rep(training_X, times = B),
      nu_train
    )

    nu_mult_Y_train <-training_Y %*% t(nu)

    real_response <- cos(as.vector(nu_mult_Y_train))
    imag_response <- sin(as.vector(nu_mult_Y_train))

    model_real <- fit_lgb(features_train,
                              real_response,
                              num_rounds,
                              params = lgb_params)
    model_imag <- fit_lgb(features_train,
                              imag_response,
                              num_rounds,
                              params = lgb_params)
    
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
    
    matrix(real_hat, nrow = N_eval, ncol = B) + 1i * matrix(imag_hat, nrow = N_eval, ncol = B)
  
}


## ----- functions to apply to time series data ------ 

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

### get pairs (t, t + 1) such that both entries are training points
pairs_from_mask <- function(in_mask) {
  Tlen <- length(in_mask)
  which(in_mask[-Tlen] & in_mask[-1])
}

#### Gets Y as a matrix
# get_Y_mat <- function(data) {
#   ycols <- grep("^Y[0-9]*$", names(data))
#   if (length(ycols) == 0 && "Y" %in% names(data)) ycols <- which(names(data) == "Y")
#   if (length(ycols) == 0) stop("No Y columns found (need 'Y' or 'Y1..Yd').")
#   as.matrix(data[, ycols, drop = FALSE])
# }
get_Y_mat <- function(data) {
  nm <- names(data)
  ycols <- grep("^Y(\\d+)?$", nm)
  if (length(ycols) == 0) stop("No Y columns found. Provide 'Y' or 'Y1..Yd'.")
  as.matrix(data[, ycols, drop = FALSE])
}