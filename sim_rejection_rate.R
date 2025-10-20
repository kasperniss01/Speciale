### document for simulating rejection rate 
### right now only works for 2D-AR(1) processes!!!

source("simulate_AR_process.R")
source("estimate_test.R")
source("sim_crit_value.R")
source("oracle_statistic.R")

sim_rej_rate <- function(Tlen, n, L, B,
                         A_matrix, 
                         alphas,
                         repetitions = 500,
                         verbose = TRUE) {
  
  mu <- matrix(rnorm(B), ncol = 1)
  nu <- matrix(rnorm(B), ncol = 1)   # d = 1 here
  K <- length(alphas)
  
  reject <- matrix(FALSE, nrow = K, ncol = repetitions,
                   dimnames = list(paste0("alpha=", alphas), NULL))
  
  reject_true <- matrix(FALSE, nrow = K, ncol = repetitions,
                   dimnames = list(paste0("alpha=", alphas), NULL))
  
  S_hats <- numeric(repetitions)
  S_trues <- numeric(repetitions)
  
  covvar_list <- list()
  covvar_true_list <- list()
  
  data_list <- list()
  
  for (i in seq_len(repetitions)) {
    data <- simulate_2D_AR1_process(Tlen, A_matrix, verbose = FALSE)
    data_list[[i]] <- data
    
    est <- estimate_stat(
      data = data, n = n, L = L, B = B,
      mu = mu, nu = nu,
      remainder_true_ccfs = list(
        true_phi = function(x, u) char_func_cond_X_next_given_X_previous_mat(A[1,1], x, u),
        true_psi = function(x, u, t) char_func_cond_Y_given_X_mat(A[1,1], A[2,1], A[2,2], t, x, u)
      )
    )
    
    # browser()
    S_hat <- est$S_hat
    S_hats[i] <- S_hat
    
    S_true <- est$S_true
    S_trues[i] <- S_true
    
    covvar_est <- est$Covvar_Est
    covvar_list[[i]] <- covvar_est
    
    covvar_est_true <- est$Covvar_Est_true
    covvar_true_list[[i]] <- covvar_est_true
    
    # browser()
    draws <- sim_crit_draws(covvar_est)
    draws_true <- sim_crit_draws(covvar_est_true)
    
    crits <- crit_from_draws(draws, alphas)
    crits_true <- crit_from_draws(draws_true, alphas)
    
    reject[, i] <- (S_hat > crits)
    reject_true[, i] <- (S_true > crits_true)
    
    #print to see how far in loop
    if (verbose && i %% 10 == 0) message(i)
  }
  
    #calculate mean and SEs
    rates <- rowMeans(reject)
    rates_true <- rowMeans(reject_true)
    
    ses   <- sqrt(rates * (1 - rates) / repetitions)
    ses_true <- sqrt(rates_true * (1 - rates_true) / repetitions)
    
    
    #store results in a dataframe with alpha, rejection rates and SEs
    rejection_rate_df <- data.frame(alpha = alphas, 
                                    rate = as.numeric(rates), 
                                    se = as.numeric(ses),
                                    rate_true = as.numeric(rates_true),
                                    se_true = as.numeric(ses_true))
    
    #store S_hat and S estimates
    estimates <- data.frame(S_hat = S_hats, S_true = S_trues)
    
    #return both dataframes
    return(list(rejection_rate_df = rejection_rate_df, 
                estimates = estimates, 
                covvar = list(est = covvar_list, true = covvar_true_list),
                data = data_list))
}




Tlen <- 1000
n <- 100
L <- 10
B <- 10

A <- matrix(c(-0.4, 0, -0.3, 0.8), nrow = 2, byrow = T)

df <- sim_rej_rate(Tlen, n, L, B, A, seq(0.01, 1, 0.01), repetitions = 10)

# df_prec <- sim_rej_rate(Tlen, n, L, B, A, seq(0.01, 1, 0.01))
# 
# df_prec_CI <- df_prec %>% mutate(lower = rate - 1.96 * se,
#                                  upper = rate + 1.96 * se)
# 
# 
# 
# ggplot(df_prec$rejection_rate_df, aes(x = alpha, y = rate)) + 
#   geom_line() + 
#   geom_abline(color = "red") + 
#   theme_bw()
# 
# ggplot(df_prec$estimates) + 
#   geom_point(aes(x = S_hat, y = S)) + 
#   geom_abline(color = "red") + 
#   theme_bw()
# 
# ggplot(df_prec$estimates %>% pivot_longer(cols = c(S_hat, S))) + 
#   geom_histogram(aes(x = value), color = "white") + 
#   facet_wrap(~name) + 
#   theme_bw()
# 
# 
# ggplot(df_prec$estimates) + 
#   geom_qq(aes(sample = S_hat),
#           distribution = function(p) quantile(df_prec$estimates$S, p)) + 
#   geom_abline(color = "red") + 
#   theme_bw()
# 
# 
# saveRDS(df_prec, file = "rej_rate_T1000.rds")



