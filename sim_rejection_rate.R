### document for simulating rejection rate 
### right now only works for 2D-AR(1) processes!!!

source("simulate_AR_process.R") #this is needed
source("estimate_test.R") #this is needed
source("sim_crit_value.R") #this is needed
source("oracle_statistic.R") #don't think this is needed
source("conditional_distributions.R") #this is needed
source("simulate_SDE.R")
source("CIR_drift_diffusion.R")
source("simulate_parameters.R")

library(vars)

#todo: clean up, new name
# make sure vars package is installed
# VERY IMPORTANT!!!!
# fix terminology and everything with AR1 vs CIR and things about A_matrix, 
# true_CCFS and parametric stuff...

#only for AR1 process, consider changing data-generating process to be argument so it works in general
sim_rej_rate <- function(Tlen, L, B,
                         DGP = c("AR1", "CIR", "IID"),
                         parameters = list(), 
                         # A_matrix, 
                         alphas,
                         repetitions = 500,
                         remainder_true_ccfs = list(true_phi = NULL, true_psi = NULL), 
                         parametric = FALSE, 
                         verbose = FALSE,
                         random_seeds = NULL) {

  # browser()
  
  DGP <- match.arg(DGP)
  
  # browser()
  # true_phi <- function(x, mu) remainder_true_ccfs$true_phi(x, mu, A)
  # true_psi <- function(x, mu, t) remainder_true_ccfs$true_psi(x, mu, A, t)
  # these are seemingly not used?
  
  if (DGP == "AR1") {
    if (is.null(parameters$A_matrix)) stop("Must provide A_matrix for AR1")
    else A_matrix = parameters$A_matrix
    
    d <- ncol(A_matrix) - 1
  }
  
  if (DGP == "CIR") {
    if (any(is.null(parameters$theta))) stop("Must provide theta1, theta2 and theta3 for CIR")
    else {
      theta <- parameters$theta
      theta1 = theta$theta1
      theta2 = theta$theta2
      theta3 = theta$theta3
    }
   
    #create CIR functions based on those parameters
    drift <- make_CIR_drift(theta1, theta2)
    diffusion <- make_CIR_diffusion(theta3)
    
    d <- ncol(theta3) - 1
    
    N <- parameters$N
    
  }
  
  if (DGP == "IID") {
    if (is.null(parameters$distribution)) stop("Must provide distribution for IID")
    else {
      distribution <- parameters$distribution
      d <- parameters$d
    }
  }
  
  # d <- ncol(A_matrix) - 1
  
  # if mu and nu should be the same for all data replications  
  # mu <- matrix(rnorm(B), ncol = 1)
  # nu <- matrix(rnorm(B * d), ncol = d)
  K <- length(alphas)
  
  
  remainders <- data.frame()
  
  S_oracle_plugins <- numeric(repetitions)
  covvar_oracle_plugin_list <- list()
  reject_oracle_plugin <- matrix(FALSE, nrow = K, ncol = repetitions,
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
  estimate_stat_obj_list <- list()
  
  
  for (i in seq_len(repetitions)) {
    if (DGP == "AR1") {
      
      
      if(!is.null(random_seeds) && length(random_seeds) == repetitions) set.seed(random_seeds[i])
      
      
      data <- simulate_AR_process(Tlen, A_matrix, d = d, verbose = verbose)
      }
    if (DGP == "CIR") {
      
      if(!is.null(random_seeds) && length(random_seeds) == repetitions) set.seed(random_seeds[i])
      
      Z0 <- rep(1, d + 1)
      time_series <- simulate_sde(Tlen = Tlen, drift = drift, diffusion = diffusion, Z0 = Z0, N = N,
                                  verbose = verbose)
      data <- time_series$discretized_path #choose path corresponding to integer values
    }
    
    if(DGP == "IID") {
      
      if(!is.null(random_seeds) && length(random_seeds) == repetitions) set.seed(random_seeds[i])
      
      data <- distribution(Tlen)
    }
    
    
    data_list[[i]] <- data
    
    # browser()
    
    # new mu and nu for each data replication
    
    if(!is.null(random_seeds) && length(random_seeds) == repetitions) set.seed(random_seeds[i])
    mu <- matrix(rnorm(B), ncol = 1)
    if(!is.null(random_seeds) && length(random_seeds) == repetitions) set.seed(random_seeds[i] + 1e4)
    nu <- matrix(rnorm(B * d), ncol = d)
    
    #change to sim_rejection_rate: if parametric, call with plugin for AR1 and true_CCFS
    # this need huge change!!!
    if (parametric) {
      est <- estimate_stat(
        data = data, L = L, B = B,
        mu = mu, nu = nu, 
        A = A_matrix, 
        parametric_plugin_AR1 = parametric,
        remainder_true_ccfs = remainder_true_ccfs
      )
    }
    else {
      est <- estimate_stat(
        data = data, L = L, B = B,
        mu = mu, nu = nu, 
        A = A_matrix 
      )
    }
    
    

    S_hat <- est$S_hat
    S_hats[i] <- S_hat
    
    covvar_est <- est$Covvar_Est
    
    all_gamma_hat <- est$all_gamma_hat
    
    
    
    covvar_list[[i]] <- covvar_est
  
    
    # draws <- sim_crit_draws(covvar_est, nrep = 1e5)
    #draws1 <- sim_crit_draws(covvar_est, nrep = 1e5)
    draws_alt <- sim_crit_draws_alt2(covvar_est, nrep = 1e5)
    
    # 
    # max(abs(crit_from_draws(sim_crit_draws(covvar_est, nrep = 1e5), alphas) - 
    #           crit_from_draws(sim_crit_draws_alt(all_gamma_hat, nrep = 1e5), alphas))) 
    # 
    # replicate(50, {
    # max(abs(crit_from_draws(sim_crit_draws_alt(all_gamma_hat, nrep = 1e5), alphas) -
    #         crit_from_draws(sim_crit_draws_alt(all_gamma_hat, nrep = 1e5), alphas)))}
    # ) -> remainders1
    # 
    # replicate(50, {
    #   max(abs(crit_from_draws(sim_crit_draws_alt2(covvar_est, nrep = 1e5), alphas) -
    #           crit_from_draws(sim_crit_draws_alt(all_gamma_hat, nrep = 1e5), alphas))) })  -> remainders2
    # 
    # tibble(method = rep(c("alt", "alt2"), each = 50),
    #        remainders = c(remainders1, remainders2)) %>%
    #   ggplot(aes(x = remainders, fill = method, color = method)) +
    #   geom_density(alpha = 0.3)

    
    
    #browser()
    
    #crits <- crit_from_draws(draws, alphas)
    #crits1 <- crit_from_draws(draws1, alphas)
    crits_alt <- crit_from_draws(draws_alt, alphas)
    
    # (crits - crits_alt %>% abs) %>% max()
    # (crits - crits1 %>% abs) %>% max()
    
    
    
    reject[, i] <- (S_hat > crits_alt)
    
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
      all_gamma_parametric_plugin <- est$parametric_plugin$all_gamma_parametric_plugin
      
      covvar_parametric_plugins_list[[i]] <- covvar_parametric_plugin
      
      #draws_parametric_plugins <- sim_crit_draws(covvar_parametric_plugin)
      
      draws_parametric_plugins <- sim_crit_draws_alt2(covvar_parametric_plugin)
      
      crits_parametric_plugins <- crit_from_draws(draws_parametric_plugins, alphas)
      
      reject_parametric_plugins[, i] <- (S_parametric_plugin > crits_parametric_plugins)
      
    }
    
    ## if calculate true stuff
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      remainders <- bind_rows(remainders, est$Remainders)
      
      S_oracle_plugin <- est$S_true
      S_oracle_plugins[i] <- S_oracle_plugin
      
      covvar_est_oracle_plugin <- est$Covvar_Est_true
      all_gamma_oracle_plugin <- est$all_gamma_true
      
      covvar_oracle_plugin_list[[i]] <- covvar_est_oracle_plugin
    
      #draws_oracle_plugin <- sim_crit_draws(covvar_est_oracle_plugin)
      draws_oracle_plugin <- sim_crit_draws_alt2(covvar_est_oracle_plugin)
      
      crits_oracle_plugin <- crit_from_draws(draws_oracle_plugin, alphas)
      reject_oracle_plugin[, i] <- (S_oracle_plugin > crits_oracle_plugin)
    }
    
    #print to see how far in loop
    #if (verbose && i %% 10 == 0) 
      message(i)
  }
  
    #calculate mean and SEs
    rates <- rowMeans(reject)
    ses   <- sqrt(rates * (1 - rates) / repetitions)
    
    
    #store results in a dataframe with alpha, rejection rates and SEs
    rejection_rate_df <- data.frame(alpha = alphas, 
                                    rate_nonparametric = as.numeric(rates), 
                                    se_nonparametric = as.numeric(ses))
    
    estimates <- data.frame(S_hat = S_hats)
    covvars <- list(est = covvar_list)
    
    if (parametric) {

      
      rates_parametric <- rowMeans(parametric_reject)
      ses_parametric <- sqrt(rates_parametric * (1 - rates_parametric) / repetitions)
      
      rates_parametric_plugin <- rowMeans(reject_parametric_plugins)
      ses_parametric_plugin <- sqrt(rates_parametric_plugin * (1 - rates_parametric_plugin) / repetitions)
      
      rejection_rate_df <- rejection_rate_df %>% mutate(
        rate_parametric = as.numeric(rates_parametric),
        se_parametric = as.numeric(ses_parametric),
        rate_parametric_plugin = as.numeric(rates_parametric_plugin),
        se_parametric_plugin = as.numeric(ses_parametric_plugin)
      )
      
      estimates <- estimates %>% mutate(S_parametric_plugin = S_parametric_plugins)
      
      covvars <- append(covvars, list(parametric_plugin = covvar_parametric_plugins_list))
    }
    
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      
      rates_oracle_plugin <- rowMeans(reject_oracle_plugin)
      ses_oracle_plugin <- sqrt(rates_oracle_plugin * (1 - rates_oracle_plugin) / repetitions)
      
      # append oracle_plugin stuff
      rejection_rate_df <- rejection_rate_df %>% mutate(
        rate_oracle_plugin = as.numeric(rates_oracle_plugin),
        se_oracle_plugin = as.numeric(ses_oracle_plugin))
      
      estimates <- estimates %>% mutate(S_oracle_plugin = S_oracle_plugins)
      
      covvars <- append(covvars, list(true = covvar_oracle_plugin_list))
    }
    
    
    output <- list(rejection_rate_df = rejection_rate_df,
                   estimates = estimates,
                   covvar = covvars,
                   data = data_list,
                   Tlen = Tlen,
                   DGP = DGP)
    
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      output <- append(output, list(remainders = remainders))
    }
    
    
    return(output)
}

### ------- testing stuff ----- ###

# d <- 4
# A_matrix <- matrix(c(0.3, 0, 0, 0, 0.2, -0.3, 0.35, 0.7, -0.4, -0.6, 0.2, 0.5, 0.2, -0.4, 0.9, 0.3), byrow = T, nrow = 4)
# theta1 <- c(0.6, 0.4, 0.4, 0.4)

# theta2 <- matrix(
#   c(1.4, 0,    0,    0,
#     0.1, 1.1, -0.05, 0,
#     0,   0.05, 0.9,  0.10,
#     0,   0,   -0.05, 0.8),
#   nrow = 4, byrow = TRUE
# )
# 
# d <- 4
# eigens <- diag(sort(runif(d, -1,1)), nrow = d)
# P_mat <- matrix(rnorm(d * d), nrow = d)
# 
# A <- solve(P_mat) %*% eigens %*% P_mat %>% round(., 2)
# 
# Anull <- A
# Anull[1, -1] <- 0
# 
# theta3 <- diag(4) * 0.5
# # 
# parameters = list(A_matrix = Anull,
#                   theta = list(theta1 = theta1, theta2 = Anull, theta3 = theta3))
# 
# Tlen <- 2000
# B <- floor(Tlen^(1/4))
# L <- 10
# 
# test_run_AR <- sim_rej_rate(Tlen, L, B, DGP = "AR1", parameters = parameters, 
#                             alphas = seq(0.01, 1, 0.01), repetitions = 200)
# 
# test_run_CIR <- sim_rej_rate(5000, 2, 2, DGP = "CIR", parameters = list(theta = theta),
#                          alphas = seq(0.01, 1, 0.01), repetitions = 200)
# 
# test_run_AR$rejection_rate_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) + 
#   geom_line() + 
#   geom_abline(color = "red") + 
#   ylim(0, 1)
# 
# test_run_CIR$rejection_rate_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) +
#   geom_line() +
#   geom_abline(color = "red") +
#   ylim(0, 1)
# 
# T10 <- 10
# L <- 2
# A <- matrix(c(0.3, 0, 0, 0.2, 0.7, -0.4, -0.6, 0.9, 0.3), byrow = T, nrow = 3)
# B <- 2
# 
# # data_test <- simulate_AR_process(Tlen, A)
# # estimate_stat(data_test, L, B)
# 
# test_df_highdim_T10 <- sim_rej_rate(10, L = 2, B = 2, DGP, A_matrix = A, alphas = seq(0.1,1, 0.05), repetitions = 10,
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
