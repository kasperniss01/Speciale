### document for simulating rejection rate 
### right now only works for 2D-AR(1) processes!!!

source("simulate_AR_process.R") #this is needed
source("estimate_stat_oracle_parametric_only.R") #this is needed
source("sim_crit_value.R") #this is needed
source("oracle_statistic.R") #don't think this is needed
source("conditional_distributions.R") #this is needed
source("simulate_SDE.R")
source("CIR_drift_diffusion.R")
source("simulate_parameters.R")

library(vars)

#todo: clean up, new name
# make sure vars package is installed
# generalize DGP

#only for AR1 process, consider changing data-generating process to be argument so it works in general
sim_rej_rate_oracle_parametric_only <- function(Tlen, L, B,
                         DGP = c("AR1", "CIR", "IID"),
                         parameters = list(), 
                         # A_matrix, 
                         alphas,
                         repetitions = 500,
                         remainder_true_ccfs = list(true_phi = NULL, true_psi = NULL), 
                         parametric = FALSE, 
                         verbose = FALSE,
                         random_seeds = seq_len(repetitions)) {
  
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
  
  # S_hats <- numeric(repetitions)
  # covvar_list <- list()
  # reject <- matrix(FALSE, nrow = K, ncol = repetitions,
  #                  dimnames = list(paste0("alpha=", alphas), NULL))
  # 
  
  S_parametric_plugins <- numeric(repetitions)
  covvar_parametric_plugins_list <- list()
  reject_parametric_plugins <- matrix(FALSE, nrow = K, ncol = repetitions,
                                      dimnames = list(paste0("alpha=", alphas), NULL))
  
  
  parametric_reject <- matrix(FALSE, nrow = K, ncol = repetitions,
                              dimnames = list(paste0("alpha=", alphas), NULL))
  
  data_list <- list()
  
  
  for (i in seq_len(repetitions)) {
    if (DGP == "AR1") {
      if(!is.null(random_seeds) && length(random_seeds) == repetitions) set.seed(random_seeds[i])
    
      data <- simulate_AR_process(Tlen, A_matrix, d = d, 
                                                  verbose = verbose)
    }
    if (DGP == "CIR") {
      Z0 <- rep(1, d + 1)
      
      if(!is.null(random_seeds) && length(random_seeds) == repetitions) set.seed(random_seeds[i])
      
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
    mu <- matrix(rnorm(B), ncol = 1)
    nu <- matrix(rnorm(B * d), ncol = d)
    
    est <- estimate_stat_oracle_parametric_only(
      data = data, L = L, B = B,
      mu = mu, nu = nu, 
      A = A_matrix, 
      parametric_plugin_AR1 = parametric,
      remainder_true_ccfs = remainder_true_ccfs
    )
    
    # S_hat <- est$S_hat
    # S_hats[i] <- S_hat
    # 
    # covvar_est <- est$Covvar_Est
    # covvar_list[[i]] <- covvar_est
    # 
    # draws <- sim_crit_draws(covvar_est)
    # crits <- crit_from_draws(draws, alphas)
    # reject[, i] <- (S_hat > crits)

    ## if parametric is TRUE
    # if (parametric) {
      # fit <- vars::VAR(data, p = 1)
      # 
      # # Extract p-values from the coefficient matrix
      # p_vals <- coefficients(fit)$X %>%
      #   t() %>%
      #   as_tibble() %>%
      #   dplyr::select(starts_with("Y")) %>%
      #   slice(4) %>%
      #   as_vector() %>%
      #   unname()
      # 
      # # Vectorized Bonferroni thresholds
      # bonferroni_alphas <- alphas / d
      # 
      # # Vectorized rejection test
      # parametric_reject[, i] <- colSums(outer(p_vals, bonferroni_alphas, "<")) > 0
      # 
      # 
      
      S_parametric_plugin <- est$parametric_plugin$S_parametric_plugin
      S_parametric_plugins[i] <- S_parametric_plugin
      
      covvar_parametric_plugin <- est$parametric_plugin$Covvar_Est_parametric_plugin
      covvar_parametric_plugins_list[[i]] <- covvar_parametric_plugin
      
      draws_parametric_plugins <- sim_crit_draws(covvar_parametric_plugin)
      crits_parametric_plugins <- crit_from_draws(draws_parametric_plugins, alphas)
      
      reject_parametric_plugins[, i] <- (S_parametric_plugin > crits_parametric_plugins)
      
      
    # }
    
    ## if calculate true stuff
    if (!is.null(remainder_true_ccfs$true_phi) && !is.null(remainder_true_ccfs$true_psi)) {
      
      S_oracle_plugin <- est$S_true
      S_oracle_plugins[i] <- S_oracle_plugin
      
      covvar_est_oracle_plugin <- est$Covvar_Est_true
      covvar_oracle_plugin_list[[i]] <- covvar_est_oracle_plugin
      
      draws_oracle_plugin <- sim_crit_draws(covvar_est_oracle_plugin)
      crits_oracle_plugin <- crit_from_draws(draws_oracle_plugin, alphas)
      reject_oracle_plugin[, i] <- (S_oracle_plugin > crits_oracle_plugin)
    }
    
    #print to see how far in loop
    #if (verbose && i %% 10 == 0) 
    message(i)
  }
  
  #calculate mean and SEs
  #rates <- rowMeans(reject)
  #ses   <- sqrt(rates * (1 - rates) / repetitions)
  
  
  #store results in a dataframe with alpha, rejection rates and SEs
  rejection_rate_df <- data.frame(alpha = alphas
                                #  rate_nonparametric = as.numeric(rates), 
                                 # se_nonparametric = as.numeric(ses)
                                )
  
  estimates <- data.frame(
   # S_hat = S_hats
    )
  
  covvars <- list(
  #  est = covvar_list
    )
  
  #if (parametric) {
    # browser()
    
    # rates_parametric <- rowMeans(parametric_reject)
    # ses_parametric <- sqrt(rates_parametric * (1 - rates_parametric) / repetitions)
    
    rates_parametric_plugin <- rowMeans(reject_parametric_plugins)
    ses_parametric_plugin <- sqrt(rates_parametric_plugin * (1 - rates_parametric_plugin) / repetitions)
    
    rejection_rate_df <- rejection_rate_df %>% mutate(
      # rate_parametric = as.numeric(rates_parametric),
      # se_parametric = as.numeric(ses_parametric),
      rate_parametric_plugin = as.numeric(rates_parametric_plugin),
      se_parametric_plugin = as.numeric(ses_parametric_plugin)
    )
    
    estimates <- data.frame(S_parametric_plugin = S_parametric_plugins)
    
    covvars <- append(covvars, list(parametric_plugin = covvar_parametric_plugins_list))
  #}
  
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
  
  
  return(output)
}




# 
# set.seed(420)
# A_mat <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)
# A_mat[1,-1] <- 0
# 
# 
# sim_rej_rate_oracle_parametric_only(
#   Tlen = 1000,
#   L = 10,
#   B = 13,
#   parameters = list(
#     A_matrix = A_mat
#   ),
#   alphas = seq(0.005, 1, 0.005),
#   remainder_true_ccfs = list(
#     true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#     true_psi = function(x, u, A, t, stationary_covariance) {stationary_ccf_of_Y_given_X(A, x, u, stationary_mean = rep(0,nrow(A)),
#                                                                                         stationary_covariance = stationary_covariance)}
#   ),
#   parametric = TRUE,
#   repetitions = 10
# )
# 






















