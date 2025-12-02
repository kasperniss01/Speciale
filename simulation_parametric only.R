source("sim_rej_rate_oracle_parametric_only.R")
source("sim_rejection_rate.R")

### simulation to see what/if rate of b-increase causes a problem for the parametric/oracle estimator.


set.seed(420)
A_mat <- runif(16, -1, 1) %>% matrix(4,4) %>% round(2)
A_mat[1,-1] <- 0



points_between_B <- 10

Tlens <- 500
Bs = c(seq(1,11, 2), seq(21,101, 20))
baseline_gamma <- c(0)

Large_rejection_rate_df2 <- data.frame()


for(B in Bs) {
  for(Tlen in Tlens) {
    
    gamma = baseline_gamma / sqrt(Tlen)
    
    Alocal <- A_mat
    Alocal[1,-1] <- gamma*rep(1/sqrt(3), 3)
    
    print(paste0("Starting Tlen = ", Tlen, ", B = ", B))
    test <- sim_rej_rate_oracle_parametric_only(
      Tlen = Tlen,
      L = 10,
      B = B,
      parameters = list(
        A_matrix = Alocal
      ),
      alphas = seq(0.005, 1, 0.005),
      remainder_true_ccfs = list(
        true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
        true_psi = function(x, u, A, t, stationary_covariance) {stationary_ccf_of_Y_given_X(A, x, u, stationary_mean = rep(0,nrow(A)), 
                                                                                            stationary_covariance = stationary_covariance)}
      ),
      parametric = TRUE,
      repetitions = 200
    )
    
    Large_rejection_rate_df2 <- rbind(Large_rejection_rate_df2,
                                     test$rejection_rate_df %>% mutate(Tlen = Tlen, B = B))
    
  }
}



Large_rejection_rate_df %>% 
  bind_rows(
    Large_rejection_rate_df1
  )  %>% 
  filter(alpha <= 0.051, alpha >= 0.049) %>% 
  ggplot(aes(x = B, y = rate_parametric_plugin)) +
  geom_line()



Large_rejection_rate_df %>% 
  bind_rows(
    Large_rejection_rate_df1
  )  %>% 
  mutate(B = factor(B)) %>% 
  ggplot(aes(x = alpha, y = rate_parametric_plugin, color = B)) + 
  geom_line() +
  geom_abline(color = "red")



#saveRDS(Large_rejection_rate_df, file = "datasets/Parametetric_plugin_B1to101_under_basegamma5.RDS")
#saveRDS(Large_rejection_rate_df, file = "datasets/Parametetric_plugin_B1to20_under_basegamma5.RDS")

Large_rejection_rate_df <- readRDS("datasets/Parametetric_plugin_B1to20_under_basegamma5.RDS")
Large_rejection_rate_df1 <- readRDS("datasets/Parametetric_plugin_B1to101_under_basegamma5.RDS")


trans_4th_root <- scales::trans_new(
  name = "4th_root",
  transform = function(x) x^(1/4),
  inverse = function(x) x^4,
  breaks = c(100,1e4))




Large_rejection_rate_df %>% 
  mutate(diff = rate_parametric_plugin - alpha, B = as.factor(B), alpha = as.factor(alpha)) %>% 
  filter(alpha %in% c("0.05", "0.2", "0.5"), B %in% c("1", "5", "10")) %>% 
  ggplot(aes(x = Tlen, y = diff, color = B)) +
  geom_line() +
  coord_trans(x = trans_4th_root) +
  scale_x_continuous(
    breaks = (3:10)^4,
    labels = (3:10)^4
  ) +
  facet_wrap(~alpha)
#  geom_vline(xintercept = (3:10)^4, linetype = "dashed", color = "grey")




Large_rejection_rate_df %>% 
  mutate(B = as.factor(B)) %>%
  filter(B == "2") %>% 
  group_by(Tlen, B) %>% 
  summarise(miscalibration = mean(abs(rate_oracle_plugin - alpha))) %>% 
  ggplot(aes(x = Tlen, y = miscalibration, color = B)) +
  geom_line() +
  coord_trans(x = trans_4th_root) +
  scale_x_continuous(
    breaks = (3:10)^4,
    labels = (3:10)^4
  )
  
  
  
  
  
  
  ggplot(aes(x = Tlen, y = rate_parametric_plugin - 0.05, color = B)) +
  geom_line() +
  coord_trans(x = trans_4th_root) +
  scale_x_continuous(
    breaks = (3:10)^4,
    labels = (3:10)^4
  )




Tlens
  
Large_rejection_rate_df %>% filter(Tlen == 80, B == 10) %>% 
  ggplot(aes(x = alpha, y = rate_parametric_plugin, color = as.factor(B))) +
  geom_line() +
  geom_abline(linetype = "dashed")



test_old <- sim_rej_rate(
  Tlen = 5000,
  L = 10,
  B = 8,
  parameters = list(
    A_matrix = A_mat
  ),
  alphas = seq(0.005, 1, 0.005),
  remainder_true_ccfs = list(
    true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
    true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
  ),
  parametric = TRUE,
  repetitions = 200
)

test_old_seed2 <- sim_rej_rate(
  Tlen = 5000,
  L = 10,
  B = 8,
  parameters = list(
    A_matrix = A_mat
  ),
  alphas = seq(0.005, 1, 0.005),
  remainder_true_ccfs = list(
    true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
    true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
  ),
  parametric = TRUE,
  repetitions = 200,
  random_seeds = seq(1001, 1001 + 200, by = 1)
)

test_old_no_fixed_seed <- sim_rej_rate(
  Tlen = 5000,
  L = 10,
  B = 8,
  parameters = list(
    A_matrix = A_mat
  ),
  alphas = seq(0.005, 1, 0.005),
  remainder_true_ccfs = list(
    true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
    true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
  ),
  parametric = TRUE,
  repetitions = 200,
  random_seeds = NULL
)





test_par_only <- sim_rej_rate_oracle_parametric_only(
  Tlen = 5000,
  L = 10,
  B = 8,
  parameters = list(
    A_matrix = A_mat
  ),
  alphas = seq(0.005, 1, 0.005),
  remainder_true_ccfs = list(
    true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
    true_psi = function(x, u, A, t, stationary_covariance) {stationary_ccf_of_Y_given_X(A, x, u, stationary_mean = rep(0,nrow(A)), 
                                                                                        stationary_covariance = stationary_covariance)}
  ),
  parametric = TRUE,
  repetitions = 200
)


test_old$rejection_rate_df %>% 
  pivot_longer(cols = c(rate_nonparametric,rate_oracle_plugin, rate_parametric_plugin), names_to = "method", values_to = "rejection_rate") %>%
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  ggtitle("Old function")

test_old_seed2$rejection_rate_df %>%
  pivot_longer(cols = c(rate_nonparametric,rate_oracle_plugin, rate_parametric_plugin), names_to = "method", values_to = "rejection_rate") %>%
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  ggtitle("Old function seed 2")

test_old_no_fixed_seed$rejection_rate_df %>% 
  pivot_longer(cols = c(rate_nonparametric,rate_oracle_plugin, rate_parametric_plugin), names_to = "method", values_to = "rejection_rate") %>%
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  ggtitle("Old function no fixed seed")


test_par_only$rejection_rate_df %>%
  pivot_longer(cols = c(rate_parametric_plugin,rate_oracle_plugin), names_to = "method", values_to = "rejection_rate") %>%
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  ggtitle("Parametric only function")



previous_run <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_0.rds")$sim_rej_obj
previous_run$rejection_rate_df %>%
  pivot_longer(cols = c(rate_nonparametric,rate_oracle_plugin, rate_parametric_plugin), names_to = "method", values_to = "rejection_rate") %>%
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_abline(linetype = "dashed") +
  ggtitle("Previous run with gamma=0")


previous_run$rejection_rate_df %>% 
  ## pivot_longer 'se_XXXX' and and 'rate_XXX' cols
  pivot_longer(cols = c(rate_nonparametric,rate_oracle_plugin, rate_parametric_plugin), 
               names_to = "method", values_to = "rejection_rate") %>%
  pivot_longer(cols = starts_with("se_"), 
               names_to = "se_method", values_to = "se_value") %>%
  # filter to only matching method and se_method
  filter((method == "rate_nonparametric" & se_method == "se_nonparametric") |
           (method == "rate_oracle_plugin" & se_method == "se_oracle_plugin") |
           (method == "rate_parametric_plugin" & se_method == "se_parametric_plugin")) %>% 
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se_value,
                  ymax = rejection_rate + 1.96 * se_value,
                  fill = method), alpha = 0.2) +
  geom_abline(linetype = "dashed") +
  ggtitle("Previous run with gamma=0 with CIs")


test_old$rejection_rate_df %>%
  ## pivot_longer 'se_XXXX' and and 'rate_XXX' cols
  pivot_longer(cols = c(rate_nonparametric,rate_oracle_plugin, rate_parametric_plugin), 
               names_to = "method", values_to = "rejection_rate") %>%
  pivot_longer(cols = starts_with("se_"), 
               names_to = "se_method", values_to = "se_value") %>%
  # filter to only matching method and se_method
  filter((method == "rate_nonparametric" & se_method == "se_nonparametric") |
           (method == "rate_oracle_plugin" & se_method == "se_oracle_plugin") |
           (method == "rate_parametric_plugin" & se_method == "se_parametric_plugin")) %>% 
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se_value,
                  ymax = rejection_rate + 1.96 * se_value,
                  fill = method), alpha = 0.2) +
  geom_abline(linetype = "dashed") +
  ggtitle("Old function with CIs")


test_par_only$rejection_rate_df %>%
  ## pivot_longer 'se_XXXX' and and 'rate_XXX' cols
  pivot_longer(cols = c(rate_parametric_plugin,rate_oracle_plugin), 
               names_to = "method", values_to = "rejection_rate") %>%
  pivot_longer(cols = starts_with("se_"), 
               names_to = "se_method", values_to = "se_value") %>%
  # filter to only matching method and se_method
  filter((method == "rate_oracle_plugin" & se_method == "se_oracle_plugin") |
           (method == "rate_parametric_plugin" & se_method == "se_parametric_plugin")) %>% 
  ggplot(aes(x = alpha, y = rejection_rate, color = method)) +
  geom_line() +
  geom_ribbon(aes(ymin = rejection_rate - 1.96 * se_value,
                  ymax = rejection_rate + 1.96 * se_value,
                  fill = method), alpha = 0.2) +
  geom_abline(linetype = "dashed") +
  ggtitle("Parametric only function with CIs")















