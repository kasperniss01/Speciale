library(tidyverse)
library(lightgbm)

simulate_AR_process <- function(Tlen, 
                                d = 1, #dimension of Y-process
                                gamma = 0, #controls H_0 or H_A
                                intercept = c(0, rep(0, d)),
                                A, #multiplied onto Z_{t-1}
                                Sigma = diag(d + 1), #variance for the noise
                                burnin = 0, #optional, to discard first samples
                                Z0 = c(0, rep(0, d)) #initial value of process
) {
  ### This functions generates samples from an AR(1) process Z = (X, Y) where 
  ### X in R and Y in R^d. The innovations are as follows
  ### Z_t = intercept + AZ_{t-1} + eps, where eps are iid noise. 
  ### The parameter gamma controls the hypothesis H_0: X_{t+1} \indep Y_t | X_t,
  ### such that gamma = 0 iff H_0 is true. 
  ### Burnin is how many samples we discard in the beginning
  
  ### Todo: implement such that innovations might be non-linear
  ### make it possible for Y to be multidimensional
  ### perform input check on dimensions for A and Sigma
  
  total_samples <- Tlen + burnin
  Z <- matrix(nrow = total_samples, ncol = d + 1) 
  Z[1, ] <- Z0
  
  #create matrix used for iid noise
  R <- chol(Sigma)
  
  #create matrix with gamma parameter - only works in Y 1D case right now
  A[1, 2] <- gamma
  
  #loop through time t
  for (t in 2:total_samples) {
    eps <- t(R) %*% rnorm(d + 1)
    Z[t, ] <- intercept + A %*% Z[t - 1, ] + eps
  }
  
  #remove burnin Z's
  Z <- Z[(burnin + 1): total_samples, ]
  
  #rename columns
  Y_name <- c()
  for(i in 1:d) Y_name[i] <- paste0("Y", i)
  colnames(Z) <- c("X", Y_name)
  as.data.frame(Z)
}

simulate_2D_AR1_process <- function(Tlen,
                                    #defaults to null-hypothesis and random numbers
                                    A_matrix = matrix(c(0.5, 0, 0.5, 0.5), nrow = 2, byrow = T),
                                    intercept = c(0, 0),
                                    Sigma = diag(2),
                                    burnin = 0,
                                    Z0 = c(0, 0),
                                    verbose = TRUE
                                    ) {
  ### This functions generates samples from an AR(1) process Z = (X, Y) in R^2 where The innovations are as follows
  ### Z_t = intercept + AZ_{t-1} + eps, where eps are iid noise. 
  
  ### A_matrix controls the null, should be lower triangular under the null
  if (verbose) print("Simulates under the hypothesis")
  
  total_samples <- Tlen + burnin
  Z <- matrix(nrow = total_samples, ncol = 2) 
  Z[1, ] <- Z0
  
  #create matrix used for iid noise
  R <- t(chol(Sigma))
  
  #loop through time t
  for (t in 2:total_samples) {
    eps <- t(R) %*% rnorm(2)
    Z[t, ] <- intercept + A_matrix %*% Z[t - 1, ] + eps
  }
  
  #remove burnin Z's
  Z <- Z[(burnin + 1): total_samples, ]
  
  #rename columns
  colnames(Z) <- c("X", "Y")
  as.data.frame(Z)
}


