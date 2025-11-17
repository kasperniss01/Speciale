### functionality to simulate parameters for AR1 and CIR processes ###

### -------------- AR1 process ----------- ###

AR1_matrix <- function(d, maxiter = 1e4, seed, gamma, 
                       e = 1/sqrt(d) * rep(1, d)) {
  #d is dimension of Y process
  #seed is for seed to generate eigenvalues
  #gamma is to control hypothesis
  #e is a unit vector
  #returns stable AR1 process matrix
  
  #simulates the eigenvalus of A to be between -1 and 1
  #creates 
  
  old_seed <- .Random.seed
  on.exit({.Random.seed <<- old_seed})
  set.seed(seed)
  
  eigens <- diag(sort(runif(d + 1, -1, 1)), nrow = d + 1)

  for(i in 1:maxiter) {
    P_mat <- matrix(rnorm((d + 1) * (d + 1)), nrow = d + 1)
    if (det(P_mat) != 0) break
  }
  
  A <- (solve(P_mat) %*% eigens %*% P_mat) %>% round(., 2)
  
  A[1, -1] <- gamma * e
  
  # browser()
  if (any((eigen(A)$values %>% abs()) >= 1)) stop("Not stable AR1 matrix")
  A
}


### ----------- CIR ----------- ###
# CIR_param <- function()

