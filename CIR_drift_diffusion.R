### functions to create lists for drift and diffusion function used as arguments to simulate_SDE
make_CIR_drift <- function(theta1, theta2) {
  #theta1 is dx1 vector
  #theta2 is dxd matrix
  
  if(all(theta2[1, -1] == 0)) {
    cat("simulating under the hypothesis \n")
    hypothesis_drift = TRUE
  }
  
  drift_fun <- function(z, t) as.numeric(theta1 - theta2 %*% z)
  
  out <- list(drift = drift_fun, hypothesis_drift = hypothesis_drift)
  
  return(out)
}

make_CIR_diffusion <- function(theta3) {
  #theta3 is dxd matrix
  #returns dxd matrix matrix diag(sqrt(z)) %*% theta3
  
  d <- nrow(theta3)
  
  if(all(theta3[1, -1] == 0)) {
    cat("simulating under the hypothesis \n")
    hypothesis_diffusion = TRUE
  }
  
  diffusion_fun <- function(z, t) diag(sqrt(z), nrow = d) %*% theta3
  
  out <- list(diffusion = diffusion_fun, hypothesis_diffusion = hypothesis_diffusion)
  
  return(out)
}