rm(list = ls())
source("sim_rejection_rate.R")

#### Run algorithm on 4D AR1 process, under local alternatives, including H0.

# B_func <- function(T) floor(T^(1/4))
L = 10
repetitions <- 200

set.seed(420)
A <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)


seeds <- 1:200


Tlens_pwr <- c(30, 100, 500, 1000)
baseline_gammas <- c(0, 5)
Bs <- c(6, 10, 30, 50, 60, 100)


for(Tlen in Tlens_pwr){
  for(baseline_gamma in baseline_gammas){
    
    if(Tlen == 30 & baseline_gamma == 5) next
    
    for (B in Bs) {
      
    
      if((baseline_gamma == 0 & B == 10)) {
        next
      }
      
      
      
      print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma, ", B = ", B))
      
      gamma = baseline_gamma / sqrt(Tlen)
      
      Alocal <- A
      Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
      
      
      simulate_temp <- sim_rej_rate(
        Tlen = Tlen,
        L = L,
        B = B,
        parameters = list(
          A_matrix = Alocal
        ),
        alphas = seq(0.005, 1, 0.005),
        remainder_true_ccfs = list(
          true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
          true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
        ),
        parametric = TRUE,
        repetitions = repetitions,
        random_seeds = seeds
      )
      
      
      sim_temp <- list(
        metadata = list(
          Tlen = Tlen,
          baseline_gamma = baseline_gamma,
          actual_gamma = gamma,
          A_matrix = Alocal,
          B = B,
          L = L,
          repetitions = repetitions
        ),
        sim_rej_obj = simulate_temp
      )
      
      
      saveRDS(sim_temp, 
              file = paste0("datasets/new_sims/VAR/L_", L, "_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B, ".rds"))
      

    }
  }
}







# 0.717
# 
# Tlens <- seq(10, 1000, by = 1)
# gamma0s <- seq(0, 20, by = 0.1)
# 
# stationary_vector <- logical(length(Tlens) * length(gamma0s))
# 
# for(i in seq_along(Tlens)) {
#   for(j in seq_along(gamma0s)) {
#     
#     gamma_temp <- gamma0s[j] / sqrt(Tlens[i])
#     
#     A_temp <- A
#     A_temp[1,-1] <- gamma_temp*rep(1/sqrt(3), 3)
#     
#     stationary_vector[(i-1)*length(gamma0s) + j] <- (max(abs(eigen(A_temp)$values)) < 1)
#   }
#}
# 
# tibble(
#   Tlen = rep(Tlens, each = length(gamma0s)),
#   gamma0 = rep(gamma0s, times = length(Tlens)),
#   stationary = stationary_vector
# ) %>% ggplot(aes(x = Tlen, y = gamma0, fill = stationary)) +
#   geom_tile() +
#   scale_fill_manual(values = c("white", "grey")) +
#   labs(title = "Stationarity region for local alternatives in 4D AR(1) process",
#        x = "Time series length T",
#        y = expression(gamma[0]),
#        fill = "Stationary") +
#   theme_minimal() +
#   geom_function(fun = function(x) sqrt(x)*0.717, color = "blue", size = 1) 
# 
# 0.717*sqrt(40)
# 






