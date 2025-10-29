rm(list = ls())
source("sim_rejection_rate.R")
source("Wald_test_AR1.R")

### 4-dimensional (i.e. 3d Y)
alternative_matrix <- function(gamma, d = 3, vec = rep(1, d), seed  = NULL) {
  vec <- vec/norm(vec, type = "2")
  
  if(!is.null(seed)) set.seed(seed)
  
  A <- matrix(runif((d+1)^2, -1,1), nrow = d+1, ncol = d+1)
  
  A <- (A / abs(eigen(A)$values[1])) %>% round(digits = 2)
  
  
  A[1,-1] <- gamma*vec
  
  A  
}

gammas <- seq(0,1,0.333)
powerTs <- c(2000, 3000)

# PowerAnalysisLists<- list()
# 
# for(Tlen in powerTs) {
#   PowerAnalysisLists[[paste0("T",as.character(Tlen))]] <- list()
#   
#   
#   for(gamma in gammas) {
#     A_alt <- alternative_matrix(gamma = gamma, d = 3, vec = c(1,1,1), seed = 100)
#     
#     message("T: ", Tlen, " | gamma: ", gamma)
#     df_power <- sim_rej_rate(Tlen, L = 10, B = 10, A = A_alt, alphas = seq(0.01, 1, 0.01),
#                             parametric = TRUE,
#                             remainder_true_ccfs = list(
#                               true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                               true_psi = function(x, u, A, t) {
#                                 char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#                               }
#                             )
#                           )
#     
#     PowerAnalysisLists[[paste0("T",as.character(Tlen))]][[paste0("gamma", as.character(gamma))]] <- df_power
#   }
#   
# }
# 
# saveRDS(PowerAnalysisLists, file = "datasets/KasperPowerAnalysis_AR1_4D.rds")




### The above gave rejection in every case for gamma > 0. 
### Now trying succesively smaller gammas, until a rejection at alpha = 0.05 is found


  # 
  # for(gamma in seq(0.33,0.001,-0.1)) {
  #   A_alt <- alternative_matrix(gamma = gamma, d = 3, vec = c(1,1,1), seed = 100)
  #   
  #   message("T: ", Tlen, " | gamma: ", gamma)
  #   df_power <- sim_rej_rate(2000, L = 10, B = 10, A = A_alt, alphas = 0.05,
  #                           parametric = TRUE,
  #                           remainder_true_ccfs = list(
  #                             true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
  #                             true_psi = function(x, u, A, t) {
  #                               char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
  #                             }
  #                           ),
  #                           repetitions = 1
  #   )
  #   print(paste0("gamma = ", gamma, " did not lead to rejection"))
  #   
  #   if(df_power$rejection_rate_df$rate < 1) break
  #   
  # }
  # 



####### gamma = 0.03 leads to rejetion. Conducting a power analysis on this gamma

# PowerAnalysis_AR1_4D$T2000$gamma0.03
#   
# df_power <- sim_rej_rate(2000, L = 10, B = 10, 
#               A = alternative_matrix(gamma = 0.03, d = 3, vec = c(1,1,1), seed = 100), 
#               alphas = seq(0.01, 1, 0.01),
#               parametric = TRUE,
#               remainder_true_ccfs = list(
#                 true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
#                 true_psi = function(x, u, A, t) {
#                   char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
#                 }
#               )
# )
# 
# 
# PowerAnalysis_AR1_4D$T2000$gamma0.03 <- df_power
# saveRDS(PowerAnalysis_AR1_4D, file = "datasets/PowerAnalysis_AR1_4D.rds")
# 
# 
# 
# PowerAnalysis_AR1_4D$T2000$gamma0.03$rejection_rate_df
# 
# 
# for(Tlen in names(PowerAnalysis_AR1_4D)) {
#   
#   for(gamma in names(PowerAnalysis_AR1_4D[[Tlen]])) {
# 
#     print(paste0(Tlen, gamma))
#     
#      PowerAnalysis_AR1_4D[[Tlen]][[gamma]] <- add_wald_test_to_sim_rej_rat_obj(PowerAnalysis_AR1_4D[[Tlen]][[gamma]])
#   }
#   
#   
# }
# 
# for(Tlen in names(PowerAnalysis_AR1_4D)) {
#   
#   for(gamma in names(PowerAnalysis_AR1_4D[[Tlen]])) {
#     
#     print(paste0(Tlen, gamma))
#     
#     PowerAnalysis_AR1_4D[[Tlen]][[gamma]]$rejection_rate_df <- PowerAnalysis_AR1_4D[[Tlen]][[gamma]]$rejection_rate_df %>% 
#         rename("rate_parametric_plugin" = "rates_parametric_plugin",
#                "rate_parametric" = "rates_parametric",
#                "rate_nonparametric" = "rate",
#                "se_nonparametric" = "se")
#   
#     
#     }
#   
#   
# }

#saveRDS(PowerAnalysis_AR1_4D, file = "datasets/PowerAnalysis_AR1_4D.rds")


PowerAnalysis_AR1_4D <- readRDS("datasets/PowerAnalysis_AR1_4D.rds")


# PLOTTING THE POWER ANALYSIS RESULTS ON T = 2000, gamma = 0.003

PowerAnalysis_AR1_4D$T2000$gamma0.03$rejection_rate_df %>%
  pivot_longer(cols = starts_with("rate"), names_to = "method", values_to = "rejection_rate") %>%
  pivot_longer(cols = starts_with("se"), names_to = "method_se", values_to = "se") %>%
  filter( str_remove(method, "rate") == str_remove(method_se, "se")) %>% 
  mutate(method = str_remove(method, "rate_")) %>%
  filter(method %in% c("parametric_plugin", "wald", "nonparametric")) %>% 
  dplyr::select(-method_se) %>%
  bind_rows(
    PowerAnalysis_AR1_4D$T2000$gamma0$rejection_rate_df %>% 
      pivot_longer(cols = starts_with("rate"), names_to = "method", values_to = "rejection_rate") %>%
      pivot_longer(cols = starts_with("se"), names_to = "method_se", values_to = "se") %>% 
      filter( str_remove(method, "rate") == str_remove(method_se, "se")) %>% 
      filter(method == "rate_nonparametric") %>% 
      mutate(method = "nonparametric under H0") %>% 
      dplyr::select(-method_se)
  ) %>% 
  
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  geom_ribbon(aes(ymin = rejection_rate - 1.96*se,
                  ymax = rejection_rate + 1.96*se,
                  fill = method), alpha = 0.2, color = NA) +
  theme_minimal() +
  labs(title = "Power Analysis in 4D AR(1) Model",
       subtitle = "Alternative with gamma = 0.03 in the first coordinate",
       y = "Rejection Rate",
       x = "Significance Level (alpha)",
       color = "Method",
       fill = "Method")




gammas <- seq(0.05,0.2,0.05)
powerTs <- c(2000)
PowerAnalysisLists_smallGamma<- list()
for(Tlen in powerTs) {
  PowerAnalysisLists_smallGamma[[paste0("T",as.character(Tlen))]] <- list()
  
  
  for(gamma in gammas) {
    A_alt <- alternative_matrix(gamma = gamma, d = 3, vec = c(1,1,1), seed = 100)
    
    message("T: ", Tlen, " | gamma: ", gamma)
    df_power <- sim_rej_rate(Tlen, L = 10, B = 10, A = A_alt, alphas = seq(0.01, 1, 0.01),
                            parametric = TRUE,
                            remainder_true_ccfs = list(
                              true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                              true_psi = function(x, u, A, t) {
                                char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                              }
                            ), repetitions = 1
    )
    
    PowerAnalysisLists_smallGamma[[paste0("T",as.character(Tlen))]][[paste0("gamma", as.character(gamma))]] <- df_power
  }
  
}


PowerAnalysisLists_smallGamma$T2000$gamma0.05$rejection_rate_df






