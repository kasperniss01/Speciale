rm(list = ls())
source("sim_rejection_rate.R")

#### Run algorithm on 4D AR1 process, under local alternatives, including H0.

Tlens_pwr <- c(5000) #c(1000, 2000, 3000, 4000, 5000)
#david kører c(0, 3, 7)
#kasper kører c(1, 5, 9)
baseline_gammas <- c(0) #c(1, 5, 9) # c(0, 1, 3, 5, 7, 9)
B_func <- function(T) floor(T^(1/4))
L = 10
repetitions <- 200

set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    
    if(Tlen == 1000 & baseline_gamma == 0){
      #skip, already done
      next
    }
    
    print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
    
    gamma = baseline_gamma / sqrt(Tlen)
    
    Alocal <- A
    Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      # B = B_func(Tlen),
      B = 7,
      parameters = list(
        A_matrix = Alocal
      ),
      alphas = seq(0.005, 1, 0.005),
      remainder_true_ccfs = list(
        true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
        true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
      ),
      parametric = TRUE,
      repetitions = repetitions
    )
    
    
    sim_temp <- list(
      metadata = list(
        Tlen = Tlen,
        baseline_gamma = baseline_gamma,
        actual_gamma = gamma,
        A_matrix = Alocal,
        B = B_func(Tlen),
        L = L,
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    saveRDS(sim_temp, file = paste0("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, ".rds"))
    
 
  }
}


simulate_temp$rejection_rate_df %>% ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  ggtitle(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma)) +
  ylab("Rejection Rate") +
  xlab("Significance Level") +
  theme_minimal()









Tlens_pwr <- c(1000, 2000, 3000, 4000, 5000)
baseline_gammas <- c(11,13,15,17)

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    
    if(Tlen == 1000 & baseline_gamma == 0){
      #skip, already done
      next
    }
    
    print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
    
    gamma = baseline_gamma / sqrt(Tlen)
    
    Alocal <- A
    Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B_func(Tlen),
      parameters = list(
        A_matrix = Alocal
      ),
      alphas = seq(0.005, 1, 0.005),
      remainder_true_ccfs = list(
        true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
        true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
      ),
      parametric = TRUE,
      repetitions = repetitions
    )
    
    
    sim_temp <- list(
      metadata = list(
        Tlen = Tlen,
        baseline_gamma = baseline_gamma,
        actual_gamma = gamma,
        A_matrix = Alocal,
        B = B_func(Tlen),
        L = L,
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    saveRDS(sim_temp, file = paste0("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, ".rds"))
    
    
  }
}


Tlens_pwr <- c(1e4)
baseline_gammas <- c(0,1,3,5,7,9,11,13,15,17)

for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    
    if(Tlen == 1000 & baseline_gamma == 0){
      #skip, already done
      next
    }
    
    print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
    
    gamma = baseline_gamma / sqrt(Tlen)
    
    Alocal <- A
    Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B_func(Tlen),
      parameters = list(
        A_matrix = Alocal
      ),
      alphas = seq(0.005, 1, 0.005),
      remainder_true_ccfs = list(
        true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
        true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
      ),
      parametric = TRUE,
      repetitions = repetitions
    )
    
    
    sim_temp <- list(
      metadata = list(
        Tlen = Tlen,
        baseline_gamma = baseline_gamma,
        actual_gamma = gamma,
        A_matrix = Alocal,
        B = B_func(Tlen),
        L = L,
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    saveRDS(sim_temp, file = paste0("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, ".rds"))
    
    
  }
}
















