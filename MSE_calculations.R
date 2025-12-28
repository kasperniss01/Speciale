set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4, 4) %>% round(2)

d <- ncol(A) - 1
B <- 10
Tlens <- c(50, 100, 200, 500, 1000, 1500, 2000, 5000, 10000)
seeds <- 1:200
train_frac <- 1000 / 1200  #fraction of data to be used for training, around 83%

results <- matrix(NA_real_, nrow = length(Tlens), ncol = 2,
                  dimnames = list(Tlens, c("MSE_phi", "MSE_psi")))

for (t_idx in seq_along(Tlens)) {
  Tlen <- Tlens[t_idx]
  
  n_train <- floor(train_frac * Tlen)
  n_total <- Tlen
  
  # store results for this Tlen over all replications
  res_T <- sapply(seeds, function(s) {
    
    set.seed(s)  # fixed seed per replication
    
    mu <- matrix(rnorm(B),   ncol = 1)
    nu <- matrix(rnorm(B*d), ncol = d)
    
    data    <- simulate_AR_process(Tlen, A) %>% as.data.frame()
    data_X  <- data$X
    data_Y  <- get_Y_mat(data)
    
    ## ----- training / evaluation indices -----
    in_train <- rep(FALSE, n_total)
    index_train <- seq_len(n_train)
    in_train[index_train] <- TRUE
    
    A_hat <- vars::VAR(data[index_train, ]) %>% Acoef() %>% .[[1]]
    
    index_eval <- (n_train + 1):n_total
    n <- length(index_eval)
    
    ## ----- phi part -----
    index_train_phi         <- pairs_from_mask(in_train)
    shifted_index_train_phi <- index_train_phi + 1
    
    X_train_phi         <- data_X[index_train_phi]
    X_shifted_train_phi <- data_X[shifted_index_train_phi]
    
    index_eval_phi         <- index_eval[1:(n - 1)]
    shifted_index_eval_phi <- index_eval_phi + 1
    
    X_eval_phi         <- data_X[index_eval_phi]
    X_shifted_eval_phi <- data_X[shifted_index_eval_phi]
    
    phi_hat_mat <- phi_hat(
      X_train_phi,
      X_shifted_train_phi,
      X_eval_phi,
      mu,
      num_rounds = 300,
      lgb_params = list()
    )
    
    phi_mat <- char_func_cond_X_next_given_X_previous_mat(
      A      = A,
      x_prev = X_eval_phi,
      u      = mu
    )
    
    phi_hat_mat_par <- char_func_cond_X_next_given_X_previous_mat(
      A      = A_hat,
      x_prev = X_eval_phi,
      u      = mu
    )
    
    SE_phi      <- abs(phi_hat_mat - phi_mat)^2
    MSE_phi     <- mean(colMeans(SE_phi))
    
    ## ----- psi part -----
    X_train_psi <- data_X[index_train]
    Y_train_psi <- data_Y[index_train, , drop = FALSE]
    
    index_eval_psi <- index_eval[1:(n - 1)]
    X_eval_psi     <- data_X[index_eval_psi]
    Y_eval_psi     <- data_Y[index_eval_psi, , drop = FALSE]
    
    psi_hat_mat <- psi_hat(
      X_train_psi,
      Y_train_psi,
      X_eval_psi,
      nu,
      num_rounds = 10,
      lgb_params = list(lambda_l2 = 2)
    )
    
    psi_mat <- stationary_ccf_of_Y_given_X(
      A                     = A,
      x_t                   = X_eval_psi,
      u                     = nu,
      stationary_covariance = stationary_covariance(A)
    )
    
    SE_psi  <- abs(psi_hat_mat - psi_mat)^2
    MSE_psi <- mean(colMeans(SE_psi))
    
    c(MSE_phi, MSE_psi)  # return a length-2 vector for this seed
  })
  
  # res_T is 2 x 200: rows = (MSE_phi, MSE_psi)
  results[t_idx, ] <- rowMeans(res_T)
  
  cat("Tlen =", Tlen,
      "-> MSE_phi =", results[t_idx, "MSE_phi"],
      ", MSE_psi =", results[t_idx, "MSE_psi"], "\n")
}

results