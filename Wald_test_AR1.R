### function to calculate a Wald statistic for the nullhypothsis (in the parametric Gaussian AR1 model)
Wald_test_AR1 <- function(data, alphas) {
  # data is dataset
  # alphas is multiple significance levels
  
  #returns a list of Wald statistic, p_values and df with reject or not
  
  theta_hat <- vars::Acoef(vars::VAR(data, p = 1))[[1]][1, -1]
  
  Z <- as.matrix(data)[-1, ]
  
  var_theta_hat <- solve(t(Z) %*% Z)[-1,-1]
  
  statistic <- t(theta_hat) %*% solve(var_theta_hat) %*% theta_hat
  
  p_value <- pchisq(statistic, df = nrow(var_theta_hat), lower.tail = FALSE)

  reject <- tibble(alpha = alphas) %>% rowwise() %>%  mutate(reject = as.numeric(p_value < alpha) )
  
  return(list(statistic = statistic, p_value = p_value, reject = reject))
  
}


### function to add rejection_rates based on a Wald test
add_wald_test_to_sim_rej_rat_obj <- function(sim_rej_rate_obj) {
  # sim_rej_rat_obj is an object of the form coming from sim_rejection_rate
  
  # returns rejection_rate df from sim_rej_rat_obj
  
  #browser()
  
  alphas <- sim_rej_rate_obj$rejection_rate_df$alpha
  
  wald_rejections <- matrix(nrow = length(alphas), ncol = length(sim_rej_rate_obj$data))
  pvals <- numeric(length(sim_rej_rate_obj$data))
  
  for(i in seq_along(sim_rej_rate_obj$data)) {
    data <- sim_rej_rate_obj$data[[i]]
    
    wald_test_result <- Wald_test_AR1(data, alphas)
    
    wald_rejections[, i] <- wald_test_result$reject$reject
    pvals[i] <- wald_test_result$p_value 
  }
  
  wald_rejections <- data.frame(alpha = alphas) %>% 
    mutate(
      rate_wald = rowMeans(wald_rejections),
      se_wald = sqrt(rate_wald * (1 - rate_wald) / length(sim_rej_rate_obj$data))
    )
  
  sim_rej_rate_obj$rejection_rate_df <- sim_rej_rate_obj$rejection_rate_df %>% 
    left_join(wald_rejections, by = "alpha")
  
  sim_rej_rate_obj
}
