### functionality to simulate parameters for AR1 and CIR processes ###

runif(1)

### -------------- AR1 process ----------- ###

# We require that the eigenvalues all lie within the unit circle #
AR1_matrix <- function(d, maxiter = 1e4, seed, gamma, 
                       e = 1/sqrt(d) * rep(1, d)) {
  #d is dimension of Y process
  #seed is for seed to generate eigenvalues
  #gamma is to control hypothesis
  #e is a unit vector
  #returns stable AR1 process matrix
  
  #simulates the eigenvalus of A to be between -1 and 1
  #creates 
  
  old_seed <- .Random.seed #this can make the function break
  on.exit({.Random.seed <<- old_seed})
  set.seed(seed)
  
  eigens <- diag(sort(runif(d + 1, -1, 1)), nrow = d + 1)

  for(i in 1:maxiter) {
    P_mat <- matrix(rnorm((d + 1) * (d + 1)), nrow = d + 1)
    if (det(P_mat) != 0) break
  }
  
  A <- (solve(P_mat) %*% eigens %*% P_mat) %>% round(2)
  
  A[1, -1] <- gamma * e
  
  # browser()
  if (any((eigen(A)$values %>% abs()) >= 1)) stop("Not stable AR1 matrix")
  A
}


### ----------- CIR ----------- ###

# requires real part of eigenvalues of theta2 is negative
# requires theta1_i > sigma_i^2 = theta3_ii
CIR_param <- function(d, sigma = rep(1, d + 1), 
                      theta1_range = c(0, 2),
                      diag_range = c(0, 2),
                      diag_margin = 0.2, 
                      maxiter = 1e4, seed, gamma, 
                      e = 1/sqrt(d) * rep(1, d)) {
  
  #set seed only for function call
  old_seed <- .Random.seed #this can make the function break
  on.exit({.Random.seed <<- old_seed})
  set.seed(seed)
  
  
  theta1 <- runif(d + 1, theta1_range[1], theta1_range[2])
  
  #make something diagonally dominant??
  theta2 <- matrix(rexp((d + 1) * (d + 1)), nrow = d + 1)
  theta2[1, 2:(d + 1)] <- gamma * e
  diag(theta2) <- 0
  kappa <- runif(d + 1, diag_range[1], diag_range[2])
  M <- diag(kappa + rowSums(theta2) + diag_margin) - theta2
  
  theta2 <- -M
  
  # browser()
  
  if(any(Re(eigen(theta2)$values) >= 0)) stop("theta2 is not stable matrix")
  
  theta3 <- diag(sigma, nrow = d + 1)
  
  
  #check Feller stuff - Niels' condition
  for (i in 1: (d + 1)) {
    while (sigma[i]^2 >= theta1[i]) theta1[i] <- runif(1, theta1_range[1], theta1_range[2])
  }
  
  #stationary mean?
  mu <- -solve(theta2, theta1)
  
  list(theta1 = theta1, theta2 = theta2, theta3 = theta3)
}

# CIR_param(d = 3, seed = 1, gamma = 0)
# runif(10)
# 
# CIR_param(d = 3, seed = 1, gamma = 1)
# runif(10)
