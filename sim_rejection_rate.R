### document for simulating rejection rate 
### right now only works for 2D-AR(1) processes!!!

source("simulate_AR_process.R")
source("estimate_test.R")
source("sim_crit_value.R")
source("oracle_statistic.R")
source("Conditional_distributions.R")

#only for AR1 process, consider changing data-generating process to be input so it works in general
sim_rej_rate <- function(Tlen, L, B,
                         A_matrix, 
                         alphas,
                         repetitions = 500,
                         remainder_true_ccfs = list(true_phi = NULL, true_psi = NULL), 
                         verbose = TRUE) {
  
  true_phi <- remainder_true_ccfs$true_phi
  true_psi <- remainder_true_ccfs$true_psi
  
  d <- ncol(A_matrix) - 1
  
  mu <- matrix(rnorm(B), ncol = 1)
  nu <- matrix(rnorm(B * d), ncol = d)
  K <- length(alphas)
  
  reject <- matrix(FALSE, nrow = K, ncol = repetitions,
                   dimnames = list(paste0("alpha=", alphas), NULL))
  
  if (!is.null(true_phi) && !is.null(true_psi)) {
    reject_true <- matrix(FALSE, nrow = K, ncol = repetitions,
                          dimnames = list(paste0("alpha=", alphas), NULL))
    
    remainders <- list()
    S_trues <- numeric(repetitions)
    
    covvar_true_list <- list()
  }
  
  S_hats <- numeric(repetitions)
  covvar_list <- list()
  data_list <- list()
  
  for (i in seq_len(repetitions)) {
    data <- simulate_AR_process(Tlen, A_matrix, d = d)
    data_list[[i]] <- data
    
    est <- estimate_stat(
      data = data, L = L, B = B,
      mu = mu, nu = nu,
      remainder_true_ccfs = remainder_true_ccfs
    )

    S_hat <- est$S_hat
    S_hats[i] <- S_hat
    
    covvar_est <- est$Covvar_Est
    covvar_list[[i]] <- covvar_est
    
    draws <- sim_crit_draws(covvar_est)
    crits <- crit_from_draws(draws, alphas)
    reject[, i] <- (S_hat > crits)
    
    ## if calculate true stuff
    if (!is.null(true_phi) && !is.null(true_psi)) {
      remainders[[i]] <- est$Remainders
      
      S_true <- est$S_true
      S_trues[i] <- S_true
      
      covvar_est_true <- est$Covvar_Est_true
      covvar_true_list[[i]] <- covvar_est_true
    
      draws_true <- sim_crit_draws(covvar_est_true)
      crits_true <- crit_from_draws(draws_true, alphas)
      reject_true[, i] <- (S_true > crits_true)
    }
    
    #print to see how far in loop
    if (verbose && i %% 10 == 0) message(i)
  }
  
    #calculate mean and SEs
    rates <- rowMeans(reject)
    ses   <- sqrt(rates * (1 - rates) / repetitions)
    
    
    #store results in a dataframe with alpha, rejection rates and SEs
    rejection_rate_df <- data.frame(alpha = alphas, 
                                    rate = as.numeric(rates), 
                                    se = as.numeric(ses))
    
    estimates <- data.frame(S_hat = S_hats)
    covvars <- list(est = covvar_list)
    
    
    if (!is.null(true_phi) && !is.null(true_psi)) {
      
      rates_true <- rowMeans(reject_true)
      ses_true <- sqrt(rates_true * (1 - rates_true) / repetitions)
      
      # append true stuff
      rejection_rate_df <- rejection_rate_df %>% mutate(
        rate_true = as.numeric(rates_true),
        se_true = as.numeric(ses_true))
      
      estimates <- estimates %>% mutate(S_true = S_trues)
      
      covvars <- append(covvars, list(true = covvar_true_list))
    }
    
    
    output <- list(rejection_rate_df = rejection_rate_df,
                   estimates = estimates,
                   covvar = covvars,
                   data = data_list)
    
    if (!is.null(true_phi) && !is.null(true_psi)) {
      output <- append(output, list(remainders = remainders))
    }
    
    return(output)
}



Tlen <- 10
L <- 2
A <- matrix(c(0.3, 0, 0, 0.2, 0.7, -0.4, -0.6, 0.9, 0.3), byrow = T, nrow = 3)
B <- 2

# data_test <- simulate_AR_process(Tlen, A)
# estimate_stat(data_test, L, B)

test_df_highdim <- sim_rej_rate(Tlen, L, B, A, c(0.1, 0.2), repetitions = 10,
                        remainder_true_ccfs = list(
                                          true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A[1,1], x, u),
                                          true_psi = function(x, u, t) {
                                            char_func_cond_Y_given_X_highdim_mat(A[1, 1], A[-1, 1], A[-1, -1], t, x, u)
                                          }
                                        ))

A_small <- matrix(c(0.3, 0, -0.2, 0.7), byrow = T, nrow = 2)
test_df_1D <- sim_rej_rate(Tlen, L, B, A_small, c(0.1, 0.2), repetitions = 10,
                           remainder_true_ccfs = list(
                             true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A[1,1], x, u),
                             true_psi = function(x, u, t) {
                               char_func_cond_Y_given_X_highdim_mat(A[1, 1], A[-1, 1], A[-1, -1], t, x, u)
                             }
                           ))



# true_psi(10, c(1,1), 1) %>% drop()
# 
# char_func_conditional_Y_given_X_highdim(A[1,1], A[-1,1], A[-1,-1], 10, 1, c(0.5, 0.4)) %>% drop





