set.seed(123)
# dXt = (6-3*Xt)*dt + 2*sqrt(Xt)*dWt
d_x <- expression( 6-3*x )
s_x <- expression( 2*sqrt(x) )

sde::sde.sim(X0=2,drift=d_x, sigma=s_x, T = 10, N = 1000) -> X
plot(X,main="Cox-Ingersoll-Ross")

d_z <- expression(6-3*z)
s_z <- expression(2*sqrt(z))

simulate_sde <- function(drift, diffusion, Z0, Tlen, N) {
  #drift and diffusion should be function of z and t - outputs vector of length(d + 1)
  #Z0 is initial value, vector of length d + 1
  
  ##todo: explain how drift and diffusion should be chosen
  
  drift_fun <- function(z, t) eval(drift)
  diffusion_fun <- function(z, t) eval(diffusion)
  
  d <- length(Z0) #dimension of the process
  Y_dim <- d - 1
  
  delta_t <- Tlen / N
  normals <- matrix(rnorm(N * d, 0, sqrt(delta_t)), nrow = N, ncol = d)
  times <- seq(0, Tlen, length.out = N + 1)
  
  out <- matrix(nrow = N + 1, ncol = d)
  out[1, ] <- Z0
  
    # browser()
  for(i in 1:N) {
    z <- out[i, ] #previous Z
    t <- times[i] #previous time
    
    drift_vec <- drift_fun(z, t)
    diffusion_val <- diffusion_fun(z, t)
    
    if (is.matrix(diffusion_val)) diffusion_mat <- diffusion_val
    else diffusion_mat <- diag(diffusion_val)
    
    out[i + 1, ] <- pmax(out[i, ] + 
      drift_vec * delta_t + as.vector(diffusion_mat %*% normals[i, ]), 0) #make sure to get positive values 
    
  }
  
  Y_name <- c()
  for(i in 1:Y_dim) Y_name[i] <- paste0("Y", i)
  names <- c("X", Y_name)
  
  return(ts(out, start = 0, deltat = delta_t, names = names))
}

my_X <- simulate_sde(d_z, s_z, c(2,2, 10, 0.5, 1), 10, 1000)
plot(my_X)
