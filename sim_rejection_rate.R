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
  
  for (i in seq_len(repetitions)) {
    data <- simulate_2D_AR1_process(Tlen, A_matrix, verbose = FALSE)
    
    est <- estimate_stat(
      data = data, n = n, L = L, B = B,
      mu = mu, nu = nu
    )
    
    S_hat <- est$S_hat
    covvar_est <- est$Covvar_Est
    
    # browser()
    draws <- sim_crit_draws(covvar_est)
    crits <- crit_from_draws(draws, alphas)
    
    reject[, i] <- (S_hat > crits)
    
    #print to see how far in loop
    if (verbose && i %% 10 == 0) message(i)
  }
  
    #calculate mean and SEs
    rates <- rowMeans(reject)
    ses   <- sqrt(rates * (1 - rates) / repetitions)
    
    #store results in a dataframe with alpha, rejection rates and SEs
    results <- data.frame(alpha = alphas, rate = as.numeric(rates), se = as.numeric(ses))
    
    return(results)
}




Tlen <- 1000
n <- 200
L <- 5
B <- 10

A <- matrix(c(-0.4, 0, -0.3, 0.8), nrow = 2, byrow = T)

df <- sim_rej_rate(Tlen, n, L, B, A, seq(0.05, 1, 0.05))

df_prec <- sim_rej_rate(Tlen, n, L, B, A, seq(0.05, 1, 0.05))




ggplot(df_prec) + 
  geom_line(aes(x = alpha, y = rate)) + 
  geom_abline(color = "red") + 
  theme_bw()


