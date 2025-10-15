simulate_AR_process <- function(n, 
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
  
  total_samples <- n + burnin
  Z <- matrix(nrow = total_samples, ncol = 2) 
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
  colnames(Z) <- c("X", "Y")
  as.data.frame(Z)
}

sim <- simulate_AR_process(gamma = 0.5, 
                           n = 1000, 
                           burnin = 50, 
                           A = matrix(c(0.1, 0.7, -0.1, 0.3), nrow = 2))

### Quick check: partial correlation ρ(X_{t+1}, Y_t | X_t) ~ 0 
rx <- lm(sim$X[2:nrow(sim)] ~ sim$X[1:(nrow(sim) - 1)])$residuals
ry <- lm(sim$Y[1:(nrow(sim) - 1)] ~ sim$X[1:(nrow(sim) - 1)])$residuals
cor(rx, ry) #near 0 under H_0?, otherwise not

### n and L so that T = nL
n <- 10
L <- 100
Time <- n + L

### function to estimate gamma
gamma_hat <- function(B, data, ...) {
  munu <- matrix(rnorm(B * 2), ncol = 2)
  mu <- munu[, 1]
  nu <- munu[, 2]
  
  X <- data$X
  Y <- data$Y
  
  I <- complex(real = 0, imaginary = 1)
  
  Gamma_hat <- list()
  
  #loop over bs
  for (b in 1:B) { 
    Gamma_hat[b] <- exp(I * nu[b] * X)
  }
  
  Gamma_hat
}






