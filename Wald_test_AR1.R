### function to calculate the Wald statistic for the nulhypothsis (in the paramatric gaussian AR1 model)
Wald_test_AR1 <- function(data, alphas) {
  
  theta_hat <- vars::Acoef(vars::VAR(data, p = 1))[[1]][1, -1]
  
  Z <- as.matrix(data)[-1, ]
  
  var_theta_hat <- solve(t(Z) %*% Z)[-1,-1]
  
  statistic <- t(theta_hat) %*% solve(var_theta_hat) %*% theta_hat
  
  p_value <- pchisq(statistic, df = nrow(var_theta_hat), lower.tail = FALSE)

  reject <- tibble(alpha = alphas) %>% rowwise() %>%  mutate(reject = as.numeric(p_value < alpha) )
  
  return(list(statistic = statistic, p_value = p_value, reject = reject))
  
}



add_wald_test_to_sim_rej_rat_obj <- function(sim_rej_rat_obj) {
  
  alphas <- sim_rej_rate_obj$rejection_rate_df$alpha
  
  wald_rejections <- matrix(nrow = length(alphas), ncol = length(sim_rej_rate_obj$data))
  pvals <- numeric(length(sim_rej_rate_obj$data))
  
  for(i in seq_along(sim_rej_rate_obj$data)) {
    #browser()
    
    data <- sim_rej_rate_obj$data[[i]]
    
    wald_test_result <- Wald_test_AR1(data, alphas)
    
    wald_rejections[, i] <- wald_test_result$reject$reject
    pvals[i] <- wald_test_result$p_value
    
    
  }
  
  #browser()
  
  wald_rejections <- data.frame(alpha = alphas) %>% 
    mutate(
      rate_wald = rowMeans(wald_rejections),
      se_wald = sqrt(rate_wald * (1 - rate_wald) / length(sim_rej_rate_obj$data))
    )
  
  
  sim_rej_rate_obj$rejection_rate_df <- sim_rej_rate_obj$rejection_rate_df %>% 
    left_join(wald_rejections, by = "alpha")
  
  sim_rej_rate_obj
}


# df_T500_4D <- add_wald_test_to_sim_rej_rat_obj(df_T500_4D)
# df_T500_4D$rejection_rate_df <- df_T500_4D$rejection_rate_df %>% rename(
#   rate_parametric = rates_parametric,
#   rate_parametric_plugin = rates_parametric_plugin,
# )
# saveRDS(df_T500_4D, file = "datasets/T500_4D.rds")
# 
# 
# df_T2K_4D <- add_wald_test_to_sim_rej_rat_obj(df_T2K_4D)
# df_T2K_4D$rejection_rate_df <- df_T2K_4D$rejection_rate_df %>% rename(
#   rate_parametric = rates_parametric,
#   rate_parametric_plugin = rates_parametric_plugin,
# )
# saveRDS(df_T2K_4D, file = "datasets/T2K_4D.rds")
# 
# df_T5K_4D <- add_wald_test_to_sim_rej_rat_obj(df_T5K_4D)
# df_T5K_4D$rejection_rate_df <- df_T5K_4D$rejection_rate_df %>% rename(
#   rate_parametric = rates_parametric,
#   rate_parametric_plugin = rates_parametric_plugin,
# )
# saveRDS(df_T5K_4D, file = "datasets/T5K_4D.rds")
# 
# df_T10K_4D <- add_wald_test_to_sim_rej_rat_obj(df_T10K_4D)
# df_T10K_4D$rejection_rate_df <- df_T10K_4D$rejection_rate_df %>% rename(
#   rate_parametric = rates_parametric,
#   rate_parametric_plugin = rates_parametric_plugin,
# )
# saveRDS(df_T10K_4D, file = "datasets/T10K_4D.rds")







