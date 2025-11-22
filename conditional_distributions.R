### conditional distributions in an AR(1) process where
### (X_t+1 , Y_t+1) = A(X_t, Y_t) + epsilon_t
### where X_t is R-valued, Y_t is R^d-valued and epsilon is iid N(0, sigma^2I_d)

### A is matrix decomposed as [[a, 0^T], [b, C]] where a is a scalar
### 0 is the R^d 0-vector, b is a R^d vector and C is dxd matrix

## linearity yields that the conditional distributions are Gaussian 
## with mean E(Y_t | X_t) and variance V(Y_t | X_t)


### TODO: consider if we want to pass noise variance
### TODO: fix function names
### TODO: above includes highdim suffix
### TODO: make sure used packages are installed


### --- WARNING: all these functions only work under the hypothesis --- ###

#install.packages("expm") #for matrix exponentiation should this be installed or?
#install netcontrol package
library(expm)

### ------------ Variance for X process --------------- ###

variance_Xt <- function(A, t, sigma1_sq = 1) {
  # A a (d+1) x (d+1) autoregressive matrix
  # t is a vector of time-inputs
  
  #returns a vector of length(t)
  
  a <- as.numeric(A[1, 1]) #a scalar 
  
  sigma1_sq * (1 - a^(2 * t)) / (1 - a^2) 
} 


### ---------- Covariance between X_t and Y_t ----------- ###

Covariance_Xt_Yt_highdim <- function(A, t, sigma1_sq = 1) {
  # A a (d+1) x (d+1) autoregressive matrix
  # t is vector of time-inputs
  
  #returns a d x length(t) matrix
  
  #TODO: fix sigma_1_sq
  
  a <- A[1, 1] #a scalar
  b <- as.vector(A[-1, 1]) #dx1 vector
  C <- as.matrix(A[-1, -1]) #dxd matrix
  
  d <- nrow(C) 
  I_d <- diag(d)
  out <- matrix(0, nrow = d, ncol = length(t))
  b <- matrix(b, ncol = 1)
  
  
  inv_Id_minus_aC <- solve(I_d - a * C)
  inv_aId_minus_C <- solve(a * I_d - C)
  
  for (i in seq_along(t)) {
    ti <- t[i]
    term1 <- (I_d - (a * C) %^% ti) %*% inv_Id_minus_aC
    term2 <- a^ti / a * ((a^ti * I_d - C %^% ti) %*% inv_aId_minus_C)
    
    out[, i] <- as.vector(a * sigma1_sq / (1 - a^2) * (term1 - term2) %*% b)
  }
  return(out)
}

### ------------------ Variance for Y process -------------- ###

variance_Yt_highdim <- function(A, t, sigma_sq = 1) {
  # A a (d+1) x (d+1) autoregressive matrix
  # t vector of time-imputs
  
  # returns d x d x length(t) array
  
  #TODO: fix sigma_sq
  
  a <- as.numeric(A[1, 1]) #a scalar
  b <- as.vector(A[-1, 1]) #dx1 vector
  C <- as.matrix(A[- 1, -1]) #dxd matrix
  
  d <- nrow(C)
  I_d <- diag(d)
  
  out <- array(0, dim = c(d, d, length(t)))
  bbt <- b %*% t(b)
  
  for(i in seq_along(t)) {
    mat <- matrix(0, nrow = d, ncol = d)
    ti <- t[i]
    
    for(j in 1:ti) {
      mat <- mat + (C %^% (j - 1)) %*% (
        bbt * variance_Xt(A, ti - j) + 
          b %*% t(Covariance_Xt_Yt_highdim(A, ti - j)) %*% t(C) + 
          C %*% Covariance_Xt_Yt_highdim(A, ti - j) %*% t(b) +
          sigma_sq * I_d
          ) %*% (t(C) %^% (j - 1))
    }
    
    out[, , i] <- mat
  }
  
  if(length(t) == 1) out <- out[, , 1, drop = FALSE] #don't know if we need this? 
  return(out)
}



variance_Yt_closed_form_array <- function(A, t, sigma1_sq = 1) {
  # A a (d+1)x(d+1) matrix
  # t a vector of time-inputs
  
  a <- as.numeric(A[1, 1])
  b <- as.vector(A[-1, 1])
  C <- as.matrix(A[-1, -1])
  
  d <- nrow(C)
  Sigma2 <- diag(d)
  
  # Initialization
  SigmaY0 <- matrix(0, d, d)  
  r0 <- rep(0, d)             
  v0 <- 0                     
  
  # --- Build block companion M ---
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
  
  M <- rbind(
    cbind(M11, M12, M13, M14),
    cbind(M21, M22, M23, M24),
    cbind(M31, M32, M33, M34),
    cbind(M41, M42, M43, M44)
  )
  
  # initial state s0 = [vec(SigmaY0); r0; v0; 1]
  s0 <- c(as.vector(SigmaY0), r0, v0, 1)
  
  # --- Eigen decomposition for fast powers M^t ---
  eg <- eigen(M)                  # may be complex; that's fine
  V  <- eg$vectors
  lam <- eg$values
  Vinv_s0 <- solve(V, s0)         # V^{-1} s0   (one solve)
  
  # matrix of λ_i^{t_j}: (dim(M) x length(t_vec))
  LamPow <- outer(lam, t, function(l, tt) l^tt)
  
  # columns are s_t for each t: s_t = V %*% diag(lam^t) %*% V^{-1} s0
  # = V %*% ((lam^t) * (V^{-1} s0))   (elementwise *)
  ST <- V %*% sweep(LamPow, 1, Vinv_s0, `*`)   # dim = dim(M) x |t_vec|
  
  # extract vec(Var(Y_t)) = first d^2 rows, then reshape to d x d x |t|
  vec_S <- ST[seq_len(d * d), , drop = FALSE]
  # ensure real (small imaginary noise can appear numerically)
  vec_S <- Re(vec_S)
  
  out <- array(NA_real_, dim = c(d, d, length(t)))
  for (j in seq_along(t)) out[, , j] <- matrix(vec_S[, j], nrow = d, ncol = d)
  
  # if a single t was passed, keep 3D array (like your first function), no drop
  out
}




### ---------- Implementation of conditional mean function --------- ###

mean_conditional_Y_given_X_highdim <- function(A, t, x_t, sigma_sq = 1) {
  # A a (d+1) x (d+1) autoregressive matrix
  # t is vector of time inputs
  # x_t is vector of x to condition on same length as t
  
  # output matrix of size d x length(t)
  
  a <- as.numeric(A[1, 1]) #scalar
  b <- as.vector(A[-1, 1]) #dx1 vector
  C <- as.matrix(A[-1, -1]) #dxd matrix
  
  cov_Xt_Yt <- Covariance_Xt_Yt_highdim(A, t, sigma_sq) # d x length(t)
  var_Xt <- variance_Xt(A, t) # t x 1
  
  scale <- x_t / var_Xt
  
  out <- t( t(cov_Xt_Yt) * scale) 
  
  return(out)
}

### --------- Implementation of conditional variance function ------ ####

variance_conditional_Y_given_X_highdim <- function(A, t, sigma_sq = 1) {
  # A a (d+1) x (d+1) autoregressive matrix
  # t is vector of time inputs
  
  # outputs 3D array of dimension d x d x length(t)
  
  a <- as.numeric(A[1, 1]) #scalar
  b <- as.vector(A[-1, 1]) #dx1 vector
  C <- as.matrix(A[-1, -1]) #dxd matrix
  
  d <- nrow(C)
  
  var_Yt <- variance_Yt_closed_form_array(A, t) #d x d x length(t)
  var_Xt <- variance_Xt(A, t) # length(t) x 1
  cov_Xt_Yt <- Covariance_Xt_Yt_highdim(A, t) # d x length(t)
  
  out <- array(0, dim = c(d, d, length(t)))
  
  for (i in seq_along(t)) {
    r <- matrix(cov_Xt_Yt[, i], ncol = 1) #d x 1
    vx <- var_Xt[i] # scalar
    
    out[, , i] <- var_Yt[, , i] - (r %*% t(r)) / vx
  }
  
  if(length(t) == 1) out <- out[, , 1, drop = FALSE] #don't know if this is needed
  
  return(out)
}


### ------ Implementation of CCFs --------- ###

### ---------- CCF of X_{t+1} | X_t ---------- ###

char_func_cond_X_next_given_X_previous_mat <- function(A, x_prev, u) {
  #A a (d+1) x (d+1) autoregressive matrix
  #x_prev a vector
  #u a vector
  
  
  #returns a matrix of dimension length(x_prev) x length(u)
  
  a <- as.numeric(A[1, 1])
  
  N <- length(x_prev)
  B <- length(u)
  
  if(is.matrix(u)) u <- drop(u)
  
  exp(1i * outer(a * x_prev, u) - 0.5 * matrix(rep(u^2, each = N), N, B))
}

### ----------- CCF of Y_t | X_t ----------- ####

char_func_cond_Y_given_X_highdim_mat <- function(A, t, x_t, u, sigma_sq = 1) {
  # A a (d+1) x (d+1) autoregressive matrix
  # t is a vector of time input
  # x_t is a sequence of x to condition on same length as t
  # u is either a matrix of size either dx1 or B x d if mu or nu is input
  
  # returns matrix of size length(t) x length(u)
  
  #browser()
  
  
  a <- as.numeric(A[1, 1]) #scalar
  b <- as.vector(A[-1, 1]) #dx1 vector
  C <- as.matrix(A[-1, -1]) #dxd matrix
  
  B <- nrow(u)
  d <- length(b)
  
  if(nrow(u) == d) u <- t(u) #transpose if u is a dx1 matrix 
  
  mean <- mean_conditional_Y_given_X_highdim(A, t, x_t, sigma_sq) # d x length(t)
  variance <- variance_conditional_Y_given_X_highdim(A, t, sigma_sq) #d x d x length(t)
  
  out <- matrix(0, nrow = length(t), ncol = B)
  for (i in seq_along(t)) {
    mean_i <- mean[, i, drop = F] # dx1
    variance_i <- variance[, , i] #dxd
    linear_term <- as.vector(u %*% mean_i) #Bx1 #throws an error right now... non-conformable arguments
    quad_term <- rowSums((u %*% variance_i) * u) # Bx1 
    out[i, ] <- exp(1i * linear_term - 0.5 * quad_term) #char. function
  }
  
  if (length(t) == 1) out <- out[1, , drop = FALSE] #don't know if needed
  
  return(out)
}


### --------- Stationary variance in VAR(1) process -------- ###
stationary_covariance <- function(A, Sigma = diag(ncol(A))) {
  netcontrol::dlyap(t(A), Sigma)
}

### equivalent to the sum characterization ###
## N(0, sum_{i = 0}^\infty A^i Sigma (A^i)^T) ##

#new variance and covariance functions#

var_Xt <- function(A) {
  #imitates variance_Xt when t is big
  W <- stationary_covariance(A)
  
  W[1, 1]
}

variance_Yt <- function(A) {
  # browser()
  d <- nrow(A)
  W <- stationary_covariance(A)
  
  W[2:d, 2:d]
}

Covariance_Xt_Yt <- function(A) {
  d <- nrow(A)
  W <- stationary_covariance(A)
  
  matrix(W[2:d, 1], ncol = 1)
}

### --------- Stationary CCF of Y | X  process -------- ###

stationary_ccf_of_Y_given_X <- function(A, x_t, u,
                                        stationary_mean = rep(0, nrow(A)), 
                                        stationary_covariance = stationary_covariance(A)
                                        ) {
 # browser()
  
  VarX <- stationary_covariance[1, 1]
  CovXY <- stationary_covariance[-1, 1, drop = F]
  VarY <- stationary_covariance[-1, -1, drop = F]
  
  mean_cond <-  x_t %*% (t(CovXY) / VarX)
  
  var_cond <- VarY - (CovXY %*% t(CovXY)) / VarX
  
  # quadratic term (N)
  quad <- rowSums((u %*% var_cond) * u)
  
  # linear term (N × T)
  Umu <- u %*% t(mean_cond)
  
  # CF(T × N)
  t( exp(1i * Umu - 0.5 * quad) )
  
}










### ----------- deprecated functions --------- ###

### For Y one-dimensional ###

# variance_Y_t <- function(a, b, c, t) {
#   b^2 / (a - c)^2 * (
#     (a^2 - a^(2 * t)) / (1 - a^2) +
#       (c^2 - c^(2 * t)) / (1 - c^2) -
#       2 * (a * c - (a * c)^t) / (1 - a * c)
#   ) + (1 - c^(2 * t)) / (1 - c^2)
# }
# 
# 
# Covariance_Xt_Yt <- function(a, b, c, t) {
#   a * b / (1 - a^2) * ((1 - (a * c)^t) / (1 - a * c) - (a^t) / (a) * (a^t - c^t) / (a - c))
# }
# 
# mean_conditional_Y_given_X <- function(a, b, c, t, x_t) {
#   cov_Xt_Yt <- Covariance_Xt_Yt(a, b, c, t)
#   var_Xt <- variance_Xt(a, t)
#   
#   cov_Xt_Yt / var_Xt * x_t
# }
# 
# variance_conditional_Y_given_X <- function(a, b, c, t) {
#   var_Yt <- variance_Y_t(a, b, c, t)
#   cov_Xt_Yt <- Covariance_Xt_Yt(a, b, c, t)
#   var_Xt <- variance_Xt(a, t)
#   
#   var_Yt - cov_Xt_Yt^2 / var_Xt
# }

# char_func_conditional_X_next_given_X_previous <- function(A, x_prev, u) {
#   #A a (d+1) x (d+1) autoregressive matrix
#   #x_prev is a vector of xs used to condition on
#   #u is a vector of evaluation points the same length as x_prev
#   
#   #returns length(u) vector
#   
#   a <- as.numeric(A[1, 1]) #scalar
#   
#   mean_X_next_given_X_prev <- a * x_prev
#   
#   exp(1i * u * mean_X_next_given_X_prev - 0.5 * u^2)
# }


## when Y is one-dimensional

 
# char_func_conditional_Y_given_X <- function(a, b, c, t, x_t, u) {
#   mean_cond <- mean_conditional_Y_given_X(a, b, c, t, x_t)
#   var_cond <- variance_conditional_Y_given_X(a, b, c, t)
#   
#   exp(1i * u * mean_cond - 0.5 * u^2 * var_cond)
# }
# 
# #the above function but returns a matrix
# char_func_cond_Y_given_X_mat <- function(a, b, c, t, x_t, u) {
#   #produces matrix
#   
#   mean_cond <- mean_conditional_Y_given_X(a, b, c, t, x_t)
#   var_cond <- variance_conditional_Y_given_X(a, b, c, t)
#   
#   if(is.matrix(u)) u <- drop(u)
#   
#   exp(1i * outer(mean_cond,  u) - 0.5 * outer(var_cond, u^2))
# }



# d = 40
# 
# variance_Yt_closed_form(a = 0.5, b= rep(0.2, d), C = matrix(1:d^2, nrow = d, ncol = d), t = 10)


# variance_Yt_closed_form <- function(A, t,
#                                     sigma1_sq = 1) {   
#   
#   a <- as.numeric(A[1, 1])
#   b <- as.vector(A[-1, 1])
#   C <- as.matrix(A[-1, -1])
#   
#   d <- nrow(C)
#   Sigma2 <- diag(d)
#   SigmaY0 <- matrix(0, d, d)
#   r0 <- rep(0, d)
#   v0 <- 0
#   
#   # Blocks for M
#   M11 <- kronecker(C, C)                                  # d^2 x d^2
#   M12 <- kronecker(b, C) + kronecker(C, b)                # d^2 x d
#   M13 <- matrix(as.vector(b %*% t(b)), nrow = d^2, ncol = 1)
#   M14 <- matrix(as.vector(Sigma2),    nrow = d^2, ncol = 1)
#   
#   M21 <- matrix(0, d, d^2);      M22 <- a * C
#   M23 <- matrix(a * b, d, 1);    M24 <- matrix(0, d, 1)
#   
#   M31 <- matrix(0, 1, d^2);      M32 <- matrix(0, 1, d)
#   M33 <- matrix(a^2, 1, 1);      M34 <- matrix(sigma1_sq, 1, 1)
#   
#   M41 <- matrix(0, 1, d^2);      M42 <- matrix(0, 1, d)
#   M43 <- matrix(0, 1, 1);        M44 <- matrix(1, 1, 1)
#   
#   # Assemble M
#   M <- rbind(
#     cbind(M11, M12, M13, M14),
#     cbind(M21, M22, M23, M24),
#     cbind(M31, M32, M33, M34),
#     cbind(M41, M42, M43, M44)
#   )
#   
#   # Initial state s0 = [vec(Sigma0); r0; v0; 1]
#   s0 <- c(as.vector(SigmaY0), as.vector(r0), v0, 1)
#   
#   # One shot: s_t = M^t s0
#   st <- (M %^% t) %*% s0
#   
#   # Extract Var(Y_t)
#   vec_Sig <- st[seq_len(d * d)]
#   matrix(vec_Sig, nrow = d, ncol = d)
# }