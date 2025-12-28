### functionality to simulate solution to SDE of the form

# dZ_t = a(Z_t, t)dt + b(Z_t, t)dW_t

#with user specific drift, a, and diffusion, b.

#see also the 'SDE' package

###

source("CIR_drift_diffusion.R")
source("simulate_parameters.R")

simulate_sde <- function(Tlen, drift, diffusion, Z0, N = NULL, 
                         verbose = FALSE) {
  #drift and diffusion should be lists with function and hypothesis
  #Z0 is initial value, vector of length d
  #Tlen is length of chain
  #N is number of grid points 
    # higher N gives better precision, defaults to 10 * Tlen
  
  #returns a time series object with the simulated solution
  
  if (is.null(N)) N <- Tlen * 10
  if (N < Tlen) stop("N must be bigger than Tlen")
  
  hypothesis_drift <- hypothesis_diffusion <- NULL 
  
  #all sorts of checks, don't know if necessary...
  #check initial value
  if (any(Z0 < 0)) stop("use positive initialization")
  # browser()
  
  #see of drift and diffusion are lists and extract stuff
  if(is.list(drift)) {
    drift_fun <- drift$drift
    hypothesis_drift <- drift$hypothesis_drift
  } 
  else stop("use make_CIR_drift to create drift function")
  
  if(is.list(diffusion)) {
    diffusion_fun <- diffusion$diffusion
    hypothesis_diffusion <- diffusion$hypothesis_diffusion
  }
  else stop("use make_CIR_diffusion to create diffusion function")
  
  if (verbose & (is.null(hypothesis_drift) | is.null(hypothesis_diffusion))) {
    cat("can't check if hypothesis is true")
  }
  else {
    if(verbose & hypothesis_drift & hypothesis_diffusion) cat("simulating under the hypothesis")
  }

  d <- length(Z0) #dimension of the process
  Y_dim <- d - 1
  
  delta_t <- Tlen / N
  normals <- matrix(rnorm(N * d, 0, sqrt(delta_t)), nrow = N, ncol = d) #scale either here or in loop
  # times <- seq(0, Tlen, length.out = N + 1)
  times <- seq(0, Tlen, by = delta_t)
  
  out <- matrix(nrow = N, ncol = d)
  out[1, ] <- Z0
  
  for(i in 1:(N - 1)) {
    z <- pmax(out[i, ], 0) #previous Z 
    # cat("iteration: ", i, "value of z:", z)
    if (any(is.na(z))) #browser()
    t <- times[i] #previous time
    
    a <- drift_fun(z, t) #drift
    b_val <- diffusion_fun(z, t) #diffusion_value
    
    #diffusion
    if (is.matrix(b_val)) b <- b_val
    else b <- diag(b_val) 
    
    out[i + 1, ] <- pmax(
      out[i, ] + a * delta_t + as.vector(b %*% normals[i, ]), 0
      ) #make sure to get positive values 
    
  }
  
  Y_name <- c()
  for(i in 1:Y_dim) Y_name[i] <- paste0("Y", i)
  names <- c("X", Y_name)
  
  #create time series object for whole path or only integer path
  ts_obj <- ts(out, start = 0, deltat = delta_t, names = names) 
  ts_obj_integers <- ts(ts_obj[cycle(ts_obj) == 1, ], start = start(ts_obj)[1], frequency = 1)
  
  out <- list(whole_path = ts_obj, discretized_path = ts_obj_integers)
  
  return(out)
}
