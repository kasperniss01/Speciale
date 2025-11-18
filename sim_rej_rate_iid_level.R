rm(list = ls())

source("sim_rejection_rate.R")

simulate_iid_gaussian_ts <- function(Tlen) {

  d <- 3
  
  mu <- rep(0, d+1)
  Sigma <- diag(d+1)
  
  data <- MASS::mvrnorm(Tlen, mu = mu, Sigma = Sigma)
  
  colnames(data) <- c("X", paste0("Y", 1:d))

  
  return(data %>% ts)  
}


B_func <- function(T) floor(T^(1/4))
L = 10
repetitions <- 200

Tlens_iid <- c(100,250,260,620,630,1290,1300,2400,2410)
for(Tlen in Tlens_iid){

    print(paste0("Tlen = ", Tlen))
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B_func(Tlen),
      parameters = list(
        distribution = simulate_iid_gaussian_ts,
        d = 3
      ),
      DGP = "IID",
      alphas = seq(0.005, 1, 0.005),
      repetitions = repetitions
    )
    
    
    sim_temp <- list(
      metadata = list(
        Tlen = Tlen,
        B = B_func(Tlen),
        L = L,
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    saveRDS(sim_temp, file = paste0("datasets/iid_level_analysis/sim_rej_rate_IID_Gaussian_4D_Tlen_", Tlen, ".rds"))
    
    
}


readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_100_baseline_gamma_0.rds")$sim_rej_obj$rejection_rate_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() +
  geom_abline() + 
  labs(title = paste0("Rejection Rate vs Alpha (Tlen = ", readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_100_baseline_gamma_0.rds")$metadata$Tlen, ")"),
       x = "Alpha",
       y = "Rejection Rate") +
  theme_minimal()


sim_temp$sim_rej_obj$rejection_rate_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() +
  geom_abline() + 
  labs(title = paste0("Rejection Rate vs Alpha (Tlen = ", sim_temp$metadata$Tlen, ")"),
       x = "Alpha",
       y = "Rejection Rate") +
  theme_minimal()

