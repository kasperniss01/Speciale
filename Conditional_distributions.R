### conditional distributions in an AR(1) process where
### (X_t+1 , Y_t+1) = A(X_t, Y_t) + epsilon_t
### where X_t is R-valued, Y_t is R^d-valued and epsilon is iid N(0, sigma^2I_d)

## linearity yields that the conditional distributions are Gaussian 
## with mean E(Y_t | X_t) and variance V(Y_t | X_t)

library(expm)

## the formula for the variance of X_t is the same for every choice of Y_t
variance_Xt <- function(a, t, sigma1_sq = 1) {
  #a is a scalar
  sigma1_sq * (1 - a^(2 * t)) / (1 - a^2) 
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
Covariance_Xt_Yt_highdim <- function(a, b, C, t, sigma1_sq = 1) {
  # a scalar
  #b d-dimensional vector
  #C dxd matrix
  #t is vector of time-inputs
  #returns a d x length(t) matrix
  
  #fix sigma1_sq later...
  
  d <- nrow(C) 
  I_d <- diag(d)
  out <- matrix(0, nrow = d, ncol = length(t))
  b <- matrix(b, ncol = 1)
  
  for (i in seq_along(t)) {
    ti <- t[i]
    term1 <- (I_d - (a * C) %^% ti) %*% solve(I_d - a * C)
    term2 <- a^ti / a * ((a^ti * I_d - C %^% ti) %*% solve(a * I_d - C))
    
    out[, i] <- as.vector(a * sigma1_sq / (1 - a^2) * (term1 - term2) %*% b)
  }
  return(out)
}

variance_Yt_highdim <- function(a, b, C, t, sigma_sq = 1) {
  # a scalar
  # b dx1 vector
  # C dxd matrix
  # t vector of time-imputs
  
  # returns d x d x length(t) array
  
  d <- nrow(C)
  I_d <- diag(d)
  
  out <- array(0, dim = c(d, d, length(t)))
  bbt <- b %*% t(b)
  
  for(i in seq_along(t)) {
    mat <- matrix(0, nrow = d, ncol = d)
    ti <- t[i]
    
    for(j in 1:ti) {
      mat <- mat + (C %^% (j - 1)) %*% (
        bbt * variance_Xt(a, ti - j) + 
          b %*% t(Covariance_Xt_Yt_highdim(a, b, C, ti - j)) %*% t(C) + 
          C %*% Covariance_Xt_Yt_highdim(a, b, C, ti - j) %*% t(b) +
          sigma_sq * I_d
          ) %*% (t(C) %^% (j - 1))
    }
    
    out[, , i] <- mat
  }
  
  if(length(t) == 1) out <- out[, , 1, drop = FALSE] #don't know if we need this?
  return(out)
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



# d = 40
# 
# variance_Yt_closed_form(a = 0.5, b= rep(0.2, d), C = matrix(1:d^2, nrow = d, ncol = d), t = 10)


mean_conditional_Y_given_X_highdim <- function(a, b, C, t, x_t, sigma_sq = 1) {
  #a is scalar
  # b is d x 1 vector
  # C is dxd matrix
  # t is vector of time inputs
  #x_t is vector of x to condition on same length as t
  
  # output matrix of size d x length(t)
  # browser()
  
  cov_Xt_Yt <- Covariance_Xt_Yt_highdim(a, b, C, t, sigma_sq) # d x length(t)
  var_Xt <- variance_Xt(a, t) # t x 1
  
  scale <- x_t / var_Xt
  
  out <- t( t(cov_Xt_Yt) * scale) 
  #why transpose this many times? to make sure that scale multiplies correctly?
  
  return(out)
}

variance_conditional_Y_given_X_highdim <- function(a, b, C, t, sigma_sq = 1) {
  #a is a scalar
  # b is dx1 vector
  # C is dxd matrix
  # t is vector of time inputs
  #outputs 3D array of dimension d x d x length(t)
  
  d <- nrow(C)
  
  var_Yt <- variance_Yt_highdim(a, b, C, t) #dxdx length(t)
  var_Xt <- variance_Xt(a, t) # length(t) x 1
  cov_Xt_Yt <- Covariance_Xt_Yt_highdim(a, b, C, t) # d x length(t)
  
  out <- array(0, dim = c(d, d, length(t)))
  
  # browser()
  for (i in seq_along(t)) {
    r <- matrix(cov_Xt_Yt[, i], ncol = 1) #d x 1
    vx <- var_Xt[i] # scalar
    
    out[, , i] <- var_Yt[, , i] - (r %*% t(r)) / vx
  }
  
  if(length(t) == 1) out <- out[, , 1, drop = FALSE] #don't know if this is needed
  
  return(out)
}

char_func_cond_Y_given_X_highdim_mat <- function(a, b, C, t, x_t, u, sigma_sq = 1) {
  # a is a scalar
  # b is a dx1 vector
  # C is a dxd matrix
  #t is a vector of time input
  #x_t is a sequence of x to condition on same length as t
  # u is either a matrix of size either dx1 or B x d if mu or nu is input
  
  B <- nrow(u)
  d <- length(b)
  
  if(nrow(u) == d) u <- t(u) #transpose if u is a dx1 matrix - why I don't know
  
  mean <- mean_conditional_Y_given_X_highdim(a, b, C, t, x_t, sigma_sq) # d x length(t)
  variance <- variance_conditional_Y_given_X_highdim(a, b, C, t, sigma_sq) #d x d x length(t)
  #variance could also just be a matrix I think if t has length 1?
  
  out <- matrix(0, nrow = length(t), ncol = B)
  for (i in seq_along(t)) {
    # browser()
    mean_i <- mean[, i, drop = F] # dx1
    variance_i <- variance[, , i] #dxd
    linear_term <- as.vector(u %*% mean_i) #Bx1
    quad_term <- rowSums((u %*% variance_i) * u) # Bx1 - why rowSums?
    out[i, ] <- exp(1i * linear_term - 0.5 * quad_term) #char. function
  }
  
  if (length(t) == 1) out <- out[1, , drop = FALSE] #dont know if needed
  
  return(out)
}










