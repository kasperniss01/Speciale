### document for simulating rejection rate 
### right now only works for 2D-AR(1) processes!!!

source("simulate_AR_process.R") #don't think this is needed
source("estimate_test.R") #this is needed
source("sim_crit_value.R") #this is needed
source("oracle_statistic.R") #don't think this is needed
source("conditional_distributions.R") #this is needed

library(vars)

#todo: clean up, new name
# make sure vars package is installed
# generalize DGP

#only for AR1 process, consider changing data-generating process to be argument so it works in general
sim_rej_rate <- function(Tlen, L, B,
                         A_matrix, 
                         alphas,
                         repetitions = 500,
                         remainder_true_ccfs = list(true_phi = NULL, true_psi = NULL), 
                         parametric = FALSE, 
                         verbose = TRUE) {
  
  # browser()
  true_phi <- remainder_true_ccfs$true_phi
  true_psi <- remainder_true_ccfs$true_psi
  
  d <- ncol(A_matrix) - 1
  
  mu <- matrix(rnorm(B), ncol = 1)
  nu <- matrix(rnorm(B * d), ncol = d)
  K <- length(alphas)
  
  reject <- matrix(FALSE, nrow = K, ncol = repetitions,
                   dimnames = list(paste0("alpha=", alphas), NULL))
  
  if (d == 1) {
    parametric_reject <- matrix(FALSE, nrow = K, ncol = repetitions,
                                dimnames = list(paste0("alpha=", alphas), NULL))
  }
  
  if (!is.null(true_phi) && !is.null(true_psi)) {
    reject_true <- matrix(FALSE, nrow = K, ncol = repetitions,
                          dimnames = list(paste0("alpha=", alphas), NULL))
    
    remainders <- data.frame()
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
    
    ## if parametric is TRUE
    if (d == 1) {
      fit <- vars::VAR(data, p = 1)
      p_val <- coefficients(fit)$X["Y1.l1", ] %>% tail(1) #only works in Y 1D case right now
      parametric_reject[, i] <- (p_val > alphas)
      
      ## make sub-category so that if parametric and phi and psi use estimated A matrix in phi and psi
    }
    
    ## if calculate true stuff
    if (!is.null(true_phi) && !is.null(true_psi)) {
      remainders <- bind_rows(remainders, est$Remainders)
      
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
    
    if (d == 1) {
      rates_parametric <- rowMeans(parametric_reject)
      ses_parametric <- sqrt(rates_parametric * (1 - rates_parametric) / repetitions)
      
      rejection_rate_df <- rejection_rate_df %>% mutate(
        rates_parametric = as.numeric(rates_parametric),
        se_parametric = as.numeric(ses_parametric)
      )
    }
    
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
                   data = data_list,
                   Tlen = Tlen)
    
    if (!is.null(true_phi) && !is.null(true_psi)) {
      output <- append(output, list(remainders = remainders))
    }
    
    return(output)
}

### ------- testing stuff ----- ###

# T10 <- 10
# L <- 2
# A <- matrix(c(0.3, 0, 0, 0.2, 0.7, -0.4, -0.6, 0.9, 0.3), byrow = T, nrow = 3)
# B <- 2
# 
# # data_test <- simulate_AR_process(Tlen, A)
# # estimate_stat(data_test, L, B)
# 
# test_df_highdim_T10 <- sim_rej_rate(T10, L, B, A, c(0.1, 0.2), repetitions = 12,
#                         remainder_true_ccfs = list(
#                                           true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                                           true_psi = function(x, u, t) {
#                                             char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#                                           }
#                                         ))
# T20 <- 20
# test_df_highdim_T20 <- sim_rej_rate(T20, L, B, A, c(0.1, 0.2), repetitions = 12,
#                                     remainder_true_ccfs = list(
#                                       true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                                       true_psi = function(x, u, t) {
#                                         char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#                                       }
#                                     ))
# 
# T100 <- 100
# test_df_highdim_T100 <- sim_rej_rate(T100, L, B, A, c(0.1, 0.2), repetitions = 12,
#                                     remainder_true_ccfs = list(
#                                       true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                                       true_psi = function(x, u, t) {
#                                         char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#                                       }
#                                     ))
# 
A_small <- matrix(c(0.3, 0, -0.2, 0.7), byrow = T, nrow = 2)
test_df_1D <- sim_rej_rate(Tlen, L, B, A_matrix = A_small, c(0.1, 0.2), repetitions = 10,
                           remainder_true_ccfs = list(
                             true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A_small, x, u),
                             true_psi = function(x, u, t) {
                               char_func_cond_Y_given_X_highdim_mat(A_small, t, x, u)
                             }
                           ))
# 
# 
# 
# # true_psi(10, c(1,1), 1) %>% drop()
# # 
# # char_func_conditional_Y_given_X_highdim(A[1,1], A[-1,1], A[-1,-1], 10, 1, c(0.5, 0.4)) %>% drop
# 
# A_power <- matrix(c(-0.2, 0.6, 0.2, -0.7), byrow = T, nrow = 2) #under alternative
# Tlen <- 100
# L <- 5
# B <- 5
# 
# test_df_power <- sim_rej_rate(Tlen, L, B,
#                               A_small, seq(0.1, 1, 0.1),
#                               repetitions = 1)
# 
# 
# fit <- vars::VAR(test_df_power$data[[1]], p = 1, type = "const")
# vars::causality(fit, cause = "Y1")$Granger
# 
# summary(fit$varresult$X)$coefficients
# summary(fit$varresult$Y)$coefficients
# 
