### Simulations for Y highdim ###

#increasing dimension of Y #

source("sim_rejection_rate.R")
source("helper_functions.R")
source("simulate_parameters.R")
Tlen <- 500
L <- 10
repetitions <- 200
seed <- 420
seeds <- 1:repetitions
B_func <- function(T) ceiling(T^(49 / 100))
# ds <- c(8, 16, 32, 64, 128)
ds <- c(256)
gamma_func <- function(d) 10/sqrt(Tlen) * sqrt(log(d) / log(4)) #gamma increases with the dimension
#gamma(4) = what is was for local alternatives

for (d in ds) {
  gammas <- c(0, gamma_func(d)) %>% round(2)
  A <- AR1_matrix(d = d, maxiter = 1e2, max_tries = 1e5, seed = seed, gamma = gammas[2])
  
  for (gamma in gammas) {
    start_time <- Sys.time()
    
    print(paste0("d = ", d, ", gamma = ", gamma))
    Alocal <- A
    Alocal[1, -1] <- gamma * rep(1 / sqrt(d), d)
    
    simulate_temp <- sim_rej_rate(
      Tlen = Tlen,
      L = L,
      B = B_func(Tlen),
      parameters = list(
        A_matrix = Alocal
      ),
      alphas = seq(0.005, 1, 0.005),
      
      repetitions = repetitions,
      random_seeds = seeds
    )
    
    sim_temp <- list(
      metadata = list(
        Tlen = Tlen,
        actual_gamma = gamma,
        A_matrix = Alocal,
        B = B_func(Tlen),
        L = L,
        repetitions = repetitions
      ),
      sim_rej_obj = simulate_temp
    )
    
    end_time <- Sys.time()
    print(paste0("time spent: ", difftime(end_time, start_time, units = "mins") %>% round(2), "minutes"))
    
    saveRDS(sim_temp,
            file = paste0("datasets/new_sims/VAR/highdim/_d_", d, "_Tlen_", Tlen, "_B_", B_func(Tlen), "_gamma_", gamma, ".rds"))
    
  }

}

# readRDS("datasets/new_sims/VAR/highdim/_500_22_0.rds") %>% comb_rej_rate_large_obj()


# ## Y begin 10D ##
# 
# # VAR #
# 
# 
# source("sim_rejection_rate.R")
# source("helper_functions.R")
# source("simulate_parameters.R")
# 
# #structure
# d <- 10
# L <- 10
# repetitions <- 200
# A <- AR1_matrix(d = d, seed = 420, gamma = 20)
# message("searching for matrix")
# Anew <- AR1_matrix(d = d, seed = 420, gamma = 100, max_tries = 1e5)
# message("search done")
# 
# B_func <- function(T) ceiling(T^(49 / 100))
# 
# #seeds
# seeds <- 1:200
# 
# 
# 
# # Tlens_pwr <- c(100, 500, 1000, 2000)
# Tlens_pwr <- c(100)
# # baseline_gammas <- c(0, 3, 5, 10, 15)
# baseline_gammas <- c(100)
# 
# for(Tlen in Tlens_pwr){
#   for(baseline_gamma in baseline_gammas){
#     # browser()
#     
#     
#     # if(Tlen == 100 & baseline_gamma > 5) next
#     
#     # if (baseline_gamma %in% c(0, 5) & Tlen != 2000) next
#     
#     print(paste0("Tlen = ", Tlen, ", baseline gamma = ", baseline_gamma))
#     
#     gamma = baseline_gamma / sqrt(Tlen)
#     
#     Alocal <- Anew
#     Alocal[1,-1] <- gamma*rep(1/sqrt(d), d)
#     
#     #quick check if A is stable matrix
#     if (any(Mod(eigen(Alocal)$values) >= 1)) next
#     
    # simulate_temp <- sim_rej_rate(
    #   Tlen = Tlen,
    #   L = L,
    #   B = B_func(Tlen),
    #   parameters = list(
    #     A_matrix = Alocal
    #   ),
    #   alphas = seq(0.005, 1, 0.005),
    #   remainder_true_ccfs = list(
    #     true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
    #     true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
    #   ),
    #   parametric = TRUE,
    #   repetitions = repetitions,
    #   random_seeds = seeds
    # )
#     
#     
    # sim_temp <- list(
    #   metadata = list(
    #     Tlen = Tlen,
    #     baseline_gamma = baseline_gamma,
    #     actual_gamma = gamma,
    #     A_matrix = Alocal,
    #     B = B_func(Tlen),
    #     L = L,
    #     repetitions = repetitions
    #   ),
    #   sim_rej_obj = simulate_temp
    # )
#     
#     saveRDS(sim_temp, 
#             file = paste0("datasets/new_sims/VAR/10D_test_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B_func(Tlen), ".rds"))
#             # file = paste0("datasets/new_sims/VAR/_Tlen_", Tlen, "_baseline_gamma_", baseline_gamma, "_B_", B_func(Tlen), ".rds"))
#     
#     
#   }
# }
# 
# readRDS("datasets/new_sims/VAR/10D_test_Tlen_100_baseline_gamma_100_B_10.rds") %>% 
#   comb_rej_rate_large_obj() %>% 
#   pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
#                values_to = "rate",
#                names_to = "Method",
#                names_prefix = "rate_") %>% 
#   mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
#   mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
#   mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
#   ggplot(aes(x = alpha, y = rate, color = Method)) + 
#   geom_line() + 
#   geom_abline(color = "darkgrey", linetype = "dashed") + 
#   labs(
#     x = expression(paste("Significance level (", alpha, ")")),
#     y = "Rejection rate") +
#   ylim(0, 1) + 
#   facet_grid(baseline_gamma ~ Tlen + B,
#              labeller = label_bquote(
#                rows = gamma[0] == .(baseline_gamma), 
#                cols = T == .(Tlen) ~ "," ~ B == .(B)  
#              )) + 
#   theme_bw() + 
#   theme(
#     strip.background = element_rect(fill = "white", color = NA),
#     legend.position = "bottom")




