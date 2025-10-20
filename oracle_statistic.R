### oracle test statistics for 2D AR(1) process

source("helper_functions.R")
source("Conditional_distributions.R")

### function only works under the null, that is gamma coef = 0
### corresponds to lower-triangular matrix
oracle_stat_2D_AR1 <- function(data, n, L, B, 
                          A_matrix, #matrix that generates the AR(1) process
                          # simulate mu and nu 
                          mu = matrix(rnorm(B), ncol = 1),
                          nu = matrix(rnorm(B), ncol = 1)
                          ) {
  #stops if not under the null
  stopifnot(A_matrix[1, 2] == 0)
  
  X <- data$X
  Y <- data$Y #perhaps use the function that greps Y-matrix
  Ymat <- if (is.matrix(Y)) Y else cbind(Y) #should not be necessary
  
  Tlen <- nrow(data)
  
  Gamma <- complex(length.out = B)
  
  # A = [[a, 0], [b, c]] - corresponds to gamma = 0
  a <- A_matrix[1, 1]
  b <- A_matrix[2, 1]
  c <- A_matrix[2, 2]
  
  mu <- as.numeric(mu)
  nu <- as.numeric(nu)
  
  for (l in 1:(L - 1)) {
    index_eval <- II_l(Tlen, n, L, l)
    
    ### stuff for phi
    index_eval_phi <- index_eval[1:(n - 1)]
    shifted_index_eval_phi <- index_eval_phi + 1
    
    X_eval_phi <- X[index_eval_phi]
    X_next_phi <- X[shifted_index_eval_phi]
    
    N_eval_phi <- length(X_eval_phi)
    
    phi <- char_func_cond_X_next_given_X_previous_mat(a, X_eval_phi, mu)
    
    # true CCF on evaluation for X
    cc_X <- exp(1i * outer(X_next_phi, mu))
    
    ### stuff for psi
    index_eval_psi <- index_eval[1:(n - 1)]
    
    X_eval_psi <- X[index_eval_psi]
    Y_eval_psi <- Ymat[index_eval_psi, , drop = F]
    
    N_eval_psi <- length(X_eval_psi)
    
    psi <- char_func_cond_Y_given_X_mat(a, b, c, index_eval_psi, X_eval_psi, nu)
    
    # true CCF on evaluation for Y
    cc_Y <- exp(1i * (Y_eval_psi %*% t(matrix(nu, ncol = 1))))
    
    ### calculate residual terms
    resX <- cc_X - phi
    resY <- cc_Y - psi
    
    lambda <- resX * resY
    
    #calculate gamma
    Gamma <- Gamma + colSums(lambda)
  }
  
  normalizer <- (n - 1) * (L - 1)
  Gamma <- Gamma / normalizer
  S <- sqrt(normalizer) * max(abs(c(Re(Gamma), Im(Gamma))))
  
  return(list(S = S))
}
