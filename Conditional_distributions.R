## conditional distribution of Y_t given X_t is normal with mean

variance_Xt <- function(a, t) {
  (1 - a^(2 * t)) / (1 - a^2) 
}


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



### Char. func.
char_func_conditional_X_next_given_X_previous <- function(a, x_prev, u) {
  mean_X_next_given_X_prev <- a * x_prev
  
  exp(1i * u * mean_X_next_given_X_prev - 0.5 * u^2)
}

#the above function but returns a matrix
char_func_cond_X_next_given_X_previous_mat <- function(a, x_prev, u) {
  #returns a matrix
  N <- length(x_prev)
  B <- length(u)
  
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
  
  exp(1i * outer(mean_cond,  u) - 0.5 * outer(var_cond, u^2))
}










