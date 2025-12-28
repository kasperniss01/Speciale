library(tidyverse)
library(lightgbm)
source("conditional_distributions.R")


#todo: make sure packages are installed
#todo: make doclines and clean up comments

simulate_AR_process <- function(Tlen, 
                                A, #multiplied onto Z_{t-1}
                                d = ncol(A) - 1, #dimension of Y-process
                                intercept = c(0, rep(0, d)),
                                Sigma = diag(d + 1), #variance for the noise
                                Z0 = MASS::mvrnorm(n = 1, #initial value
                                                   mu = rep(0, d + 1), 
                                                   Sigma = stationary_covariance(A, Sigma)),
                                different_noiseprocess = NULL,
                                verbose = FALSE 
) {
  ### This functions generates samples from an AR(1) process Z = (X, Y) where 
  ### X in R and Y in R^(d-1). The innovations are as follows
  ### Z_t = intercept + AZ_{t-1} + eps
  
  ### perform input check on dimensions for A and Sigma
  
  if (verbose) {
    if (all(A[1, -1] == 0)) cat("simulating under the hypothesis \n")
  }
  
  
  
  Z <- matrix(nrow = Tlen, ncol = d + 1) 
  Z[1, ] <- Z0
  
  
  if(is.null(different_noiseprocess)) {
    
  #create matrix used for iid noise
  R <- chol(Sigma)
  #create matrix with gamma parameter - only works in Y 1D case right now
  # A[1, 2] <- gamma
  #loop through time t
  for (t in 2:Tlen) {
    eps <- t(R) %*% rnorm(d + 1)
    Z[t, ] <- intercept + A %*% Z[t - 1, ] + eps
  }
  } else {
    for (t in 2:Tlen) {
      eps <- different_noiseprocess(d+1)
      Z[t, ] <- intercept + A %*% Z[t - 1, ] + eps
    } 
  }
  
  #rename columns
  Y_name <- c()
  for(i in 1:d) Y_name[i] <- paste0("Y", i)
  
  # create time-series object
  out <- ts(Z, start = 0, frequency = 1, names = c("X", Y_name))
  
  return(out)
}
