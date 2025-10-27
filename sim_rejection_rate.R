### document for simulating rejection rate 
### right now only works for 2D-AR(1) processes!!!
rm(list = ls())


source("simulate_AR_process.R") #this is needed
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
  # true_phi <- function(x, mu) remainder_true_ccfs$true_phi(x, mu, A)
  # true_psi <- function(x, mu, t) remainder_true_ccfs$true_psi(x, mu, A, t)
  # these are seemingly not used?
  
  
  d <- ncol(A_matrix) - 1
  
  mu <- matrix(rnorm(B), ncol = 1)
  nu <- matrix(rnorm(B * d), ncol = d)
  K <- length(alphas)
  
  
  remainders <- data.frame()
  
  S_trues <- numeric(repetitions)
  covvar_true_list <- list()
  reject_true <- matrix(FALSE, nrow = K, ncol = repetitions,
                        dimnames = list(paste0("alpha=", alphas), NULL))
  
  S_hats <- numeric(repetitions)
  covvar_list <- list()
  reject <- matrix(FALSE, nrow = K, ncol = repetitions,
                   dimnames = list(paste0("alpha=", alphas), NULL))
  
  
  S_parametric_plugins <- numeric(repetitions)
  covvar_parametric_plugins_list <- list()
  reject_parametric_plugins <- matrix(FALSE, nrow = K, ncol = repetitions,
                                      dimnames = list(paste0("alpha=", alphas), NULL))
  
  
  parametric_reject <- matrix(FALSE, nrow = K, ncol = repetitions,
                              dimnames = list(paste0("alpha=", alphas), NULL))
  
  data_list <- list()
  
  
  for (i in seq_len(repetitions)) {
    data <- simulate_AR_process(Tlen, A_matrix, d = d)
    data_list[[i]] <- data
    
    # browser()
    
    est <- estimate_stat(
      data = data, L = L, B = B,
      mu = mu, nu = nu, 
      A = A_matrix, 
      parametric_plugin_AR1 = parametric,
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
    if (parametric) {
      fit <- vars::VAR(data, p = 1)
      
      # Extract p-values from the coefficient matrix
      p_vals <- coefficients(fit)$X %>%
        t() %>%
        as_tibble() %>%
        dplyr::select(starts_with("Y")) %>%
        slice(4) %>%
        as_vector() %>%
        unname()
      
      # Vectorized Bonferroni thresholds
      bonferroni_alphas <- alphas / d
      
      # Vectorized rejection test
      parametric_reject[, i] <- colSums(outer(p_vals, bonferroni_alphas, "<")) > 0
      
      
      #browser()
      
      S_parametric_plugin <- est$parametric_plugin$S_parametric_plugin
      S_parametric_plugins[i] <- S_parametric_plugin
      
      covvar_parametric_plugin <- est$parametric_plugin$Covvar_Est_parametric_plugin
      covvar_parametric_plugins_list[[i]] <- covvar_parametric_plugin
      
      draws_parametric_plugins <- sim_crit_draws(covvar_parametric_plugin)
      crits_parametric_plugins <- crit_from_draws(draws_parametric_plugins, alphas)
      
      reject_parametric_plugins[, i] <- (S_parametric_plugin > crits)
      
    }
    
    ## if calculate true stuff
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
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
    
    if (parametric) {
     # browser()
      
      rates_parametric <- rowMeans(parametric_reject)
      ses_parametric <- sqrt(rates_parametric * (1 - rates_parametric) / repetitions)
      
      rates_parametric_plugin <- rowMeans(reject_parametric_plugins)
      ses_parametric_plugin <- sqrt(rates_parametric_plugin * (1 - rates_parametric_plugin) / repetitions)
      
      rejection_rate_df <- rejection_rate_df %>% mutate(
        rates_parametric = as.numeric(rates_parametric),
        se_parametric = as.numeric(ses_parametric),
        rates_parametric_plugin = as.numeric(rates_parametric_plugin),
        se_parametric_plugin = as.numeric(ses_parametric_plugin)
      )
      
      estimates <- estimates %>% mutate(S_parametric_plugin = S_parametric_plugins)
      
      covvars <- append(covvars, list(parametric_plugin = covvar_parametric_plugins_list))
    }
    
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      
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
    
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      output <- append(output, list(remainders = remainders))
    }
    
    
    return(output)
}

### ------- testing stuff ----- ###

#T10 <- 10
#L <- 2
# A <- matrix(c(0.3, 0, 0, 0.2, 0.7, -0.4, -0.6, 0.9, 0.3), byrow = T, nrow = 3)
#B <- 2
# 
# # data_test <- simulate_AR_process(Tlen, A)
# # estimate_stat(data_test, L, B)
# 
# test_df_highdim_T10 <- sim_rej_rate(10, L = 2, B = 2, A_matrix = A, alphas = seq(0.1,1, 0.05), repetitions = 10,
#                         remainder_true_ccfs = list(
#                                           true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                                           true_psi = function(x, u, A, t) {
#                                             char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#                                           }
#                                         ),
#                         parametric = TRUE)


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
# A_small <- matrix(c(0.3, 0, -0.2, 0.7), byrow = T, nrow = 2)
# #test_data <- simulate_AR_process(Tlen = 500, A = A_small)
# test_df_1D <- sim_rej_rate(100, L = 10, B = 10, A_matrix = A_small, seq(0.01, 1, 0.01), repetitions = 100,
#                            parametric = T,
#                            remainder_true_ccfs = list(
#                              true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                              true_psi = function(x, u, A, t) {
#                                char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#                              }
#                            ))
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
