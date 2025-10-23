### conditional distributions in an AR(1) process where
### (X_t+1 , Y_t+1) = A(X_t, Y_t) + epsilon_t
### where X_t is R-valued, Y_t is R^d-valued and epsilon is iid N(0, sigma^2I_d)

## linearity yields that the conditional distributions are Gaussian 
## with mean E(Y_t | X_t) and variance V(Y_t | X_t)

library(expm)

## the formula for the variance of X_t is the same for every choice of Y_t
variance_Xt <- function(a, t) {
  #a should be a scalar in all cases
  (1 - a^(2 * t)) / (1 - a^2) 
}

### For Y one-dimensional
variance_Y_t <- function(a, b, c, t) {
  b^2 / (a - c)^2 * (
    (a^2 - a^(2 * t)) / (1 - a^2) +
    (c^2 - c^(2 * t)) / (1 - c^2) -
    2 * (a * c - (a * c)^t) / (1 - a * c)
  ) + (1 - c^(2 * t)) / (1 - c^2)
}


Covariance_Xt_Yt <- function(a, b, c, t) {
  a * b / (1 - a^2) * ((1 - (a * c)^t) / (1 - a * c) - (a^t) / (a) * (a^t - c^t) / (a - c))
}

mean_conditional_Y_given_X <- function(a, b, c, t, x_t) {
  cov_Xt_Yt <- Covariance_Xt_Yt(a, b, c, t)
  var_Xt <- variance_Xt(a, t)
  
  cov_Xt_Yt / var_Xt * x_t
}

variance_conditional_Y_given_X <- function(a, b, c, t) {
  var_Yt <- variance_Y_t(a, b, c, t)
  cov_Xt_Yt <- Covariance_Xt_Yt(a, b, c, t)
  var_Xt <- variance_Xt(a, t)
  
  var_Yt - cov_Xt_Yt^2 / var_Xt
}


## for d-dimensional Y
Covariance_Xt_Yt_highdim <- function(a, b, C, t, sigma_sq = 1) {
  # a scalar
  #b d-dimensional vector
  #C dxd matrix
  #returns a dx1 vector
  
  d <- nrow(C) 
  I_d <- diag(d)
  
  term1 <- (I_d - (a * C) %^% t) %*% solve(I_d - a * C)
  term2 <- a^t / a * ((a^t * I_d - C %^% t) %*% solve(a * I_d - C))
  
  a * sigma_sq / (1 - a^2) * (term1 - term2) %*% b
}

variance_Yt_highdim <- function(a, b, C, t, sigma_sq = 1) {
  # a scalar
  # b dx1 vector
  # C dxd matrix
  
  d <- nrow(C)
  I_d <- diag(d)
  mat <- matrix(0, nrow = d, ncol = d)
  
  bbt <- b %*% t(b)
  
  for(i in 1:t) {
    mat <- mat + (C %^% (i - 1)) %*% (
      bbt * variance_Xt(a, t - i) + 
        b %*% t(Covariance_Xt_Yt_highdim(a, b, C, t - i)) %*% t(C) + 
        C %*% Covariance_Xt_Yt_highdim(a, b, C, t - i) %*% t(b) +
        sigma_sq * I_d
        ) %*% (t(C) %^% (i - 1))
  }
  return(mat)
}

# Fast - chat solution
variance_Yt_closed_form <- function(a, b, C, t,
                                    sigma1_sq = 1,           # Var(eps1)
                                    Sigma2 = NULL,           # Var(eps2), default I_d
                                    SigmaY0 = NULL,          # Var(Y_0)
                                    v0 = 0,                  # Var(X_0)
                                    r0 = NULL) {             # Cov(Y_0, X_0)
  d <- nrow(C)
  if (is.null(Sigma2))  Sigma2 <- diag(d)
  if (is.null(SigmaY0)) SigmaY0 <- matrix(0, d, d)
  if (is.null(r0))      r0 <- rep(0, d)
  
  # Blocks for M
  M11 <- kronecker(C, C)                                  # d^2 x d^2
  M12 <- kronecker(b, C) + kronecker(C, b)                # d^2 x d
  M13 <- matrix(as.vector(b %*% t(b)), nrow = d^2, ncol = 1)
  M14 <- matrix(as.vector(Sigma2),    nrow = d^2, ncol = 1)
  
  M21 <- matrix(0, d, d^2);      M22 <- a * C
  M23 <- matrix(a * b, d, 1);    M24 <- matrix(0, d, 1)
  
  M31 <- matrix(0, 1, d^2);      M32 <- matrix(0, 1, d)
  M33 <- matrix(a^2, 1, 1);      M34 <- matrix(sigma1_sq, 1, 1)
  
  M41 <- matrix(0, 1, d^2);      M42 <- matrix(0, 1, d)
  M43 <- matrix(0, 1, 1);        M44 <- matrix(1, 1, 1)
  
  # Assemble M
  M <- rbind(
    cbind(M11, M12, M13, M14),
    cbind(M21, M22, M23, M24),
    cbind(M31, M32, M33, M34),
    cbind(M41, M42, M43, M44)
  )
  
  # Initial state s0 = [vec(Sigma0); r0; v0; 1]
  s0 <- c(as.vector(SigmaY0), as.vector(r0), v0, 1)
  
  # One shot: s_t = M^t s0
  st <- (M %^% t) %*% s0
  
  # Extract Var(Y_t)
  vec_Sig <- st[seq_len(d * d)]
  matrix(vec_Sig, nrow = d, ncol = d)
}




### Char. func.

## when Y is one-dimensional
char_func_conditional_X_next_given_X_previous <- function(a, x_prev, u) {
  mean_X_next_given_X_prev <- a * x_prev
  
  exp(1i * u * mean_X_next_given_X_prev - 0.5 * u^2)
}

#the above function but returns a matrix
char_func_cond_X_next_given_X_previous_mat <- function(a, x_prev, u) {
  #returns a matrix
  N <- length(x_prev)
  B <- length(u)
  
  if(is.matrix(u)) u <- drop(u)
  
  exp(1i * outer(a * x_prev, u) - 0.5 * matrix(rep(u^2, each = N), N, B))
}

char_func_conditional_Y_given_X <- function(a, b, c, t, x_t, u) {
  mean_cond <- mean_conditional_Y_given_X(a, b, c, t, x_t)
  var_cond <- variance_conditional_Y_given_X(a, b, c, t)
  
  exp(1i * u * mean_cond - 0.5 * u^2 * var_cond)
}

#the above function but returns a matrix
char_func_cond_Y_given_X_mat <- function(a, b, c, t, x_t, u) {
  #produces matrix
  
  mean_cond <- mean_conditional_Y_given_X(a, b, c, t, x_t)
  var_cond <- variance_conditional_Y_given_X(a, b, c, t)
  
  if(is.matrix(u)) u <- drop(u)
  
  exp(1i * outer(mean_cond,  u) - 0.5 * outer(var_cond, u^2))
}










