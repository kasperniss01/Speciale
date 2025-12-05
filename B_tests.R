
comb_rej_rate_large_obj <- function(...) {
  # browser()
  arguments <- list(...)
  number_dfs <- length(arguments)
  
  df <- data.frame()
  
  for (i in 1:number_dfs) {
    arg_i <- arguments[[i]]
    metadata_i <- arg_i$metadata
    df_i <- arg_i$sim_rej_obj$rejection_rate_df
    
    df_i <- df_i %>% mutate(Tlen = metadata_i$Tlen,
                            B = metadata_i$B,
                            baseline_gamma = metadata_i$baseline_gamma,
                            actual_gamma = metadata_i$actual_gamma)
    
    df <- bind_rows(df, df_i)
  }
  
  df
}


Importance_of_B_VAR1 <- comb_rej_rate_large_obj(
 readRDS("datasets/new_sims/VAR/L_10_Tlen_30_baseline_gamma_0_B_6.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_30_B_10.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_30_baseline_gamma_0_B_30.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_30_baseline_gamma_0_B_60.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_30_baseline_gamma_0_B_100.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_0_B_6.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_B_10.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_0_B_30.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_0_B_60.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_0_B_100.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_5_B_6.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_5_B_10.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_5_B_30.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_5_B_60.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_5_B_100.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_0_B_6.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_B_10.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_0_B_30.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_0_B_60.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_0_B_100.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_5_B_6.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_5_B_10.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_5_B_30.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_5_B_60.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_5_B_100.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_0_B_6.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_B_10.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_0_B_30.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_0_B_60.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_0_B_100.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_5_B_6.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_5_B_10.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_5_B_30.rds"),
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_5_B_60.rds")
 readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_5_B_100.rds")
)


Importance_of_B_VAR1 %>% 
  mutate(B = factor(B)) %>% 
  pivot_longer(cols = starts_with("rate_"), 
               names_to = "method", 
               values_to = "rate",
               names_prefix = "rate_") %>%
  filter(method == "nonparametric"
      #   ,  B %in% c(10,30,60,100)
  ) %>% 
  ggplot(aes(x = alpha, y = rate, color = B)) +
  geom_line() +
  facet_grid(
    rows = vars(baseline_gamma),
    cols = vars(Tlen),
    labeller = labeller(
      baseline_gamma = function(x) paste0("Baseline gamma: ", x),
      Tlen = function(x) paste0("Time series length: ", x)
    )
  ) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  theme_minimal() +
  coord_fixed() |
Importance_of_B_VAR1 %>% 
  mutate(lab = "Among all alpha") %>%
  bind_rows(
    Importance_of_B_VAR1 %>%
      filter(alpha <= 0.10) %>%
      mutate(lab = "Among alpha <= 0.10")
  ) %>%
  mutate(diff = rate_nonparametric - alpha) %>%
  group_by(Tlen, baseline_gamma, B, lab) %>%
  summarize(max_diff = max(diff)) %>%
  mutate(baseline_gamma = factor(baseline_gamma)) %>% 
  ggplot(aes(x = B, y = max_diff, color = baseline_gamma)) +
  geom_point() +
  geom_line() +
  facet_grid(
    cols = vars(Tlen), 
    rows = vars(lab),
    labeller = labeller(
      Tlen = function(x) paste0("Time series length: ", x)
    )) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  theme_bw() +
  ylab("Maximum difference between rejection rate and alpha")




Importance_of_B_VAR1 %>% 
  filter(alpha >= 0.049, alpha <= 0.051) %>%
  pivot_longer(cols = starts_with("rate_"), 
               names_to = "method", 
               values_to = "rate",
               names_prefix = "rate_") %>%
  filter(method == "nonparametric") %>%
  mutate(baseline_gamma = factor(baseline_gamma)) %>%
  ggplot(aes(x = B, y = rate, color = baseline_gamma)) +
  geom_line() +
  geom_point() +
  geom_abline(slope = 0, intercept = 0.05, linetype = "dashed", color = "grey") +
  facet_grid(
    cols = vars(Tlen),
    labeller = labeller(
      baseline_gamma = function(x) paste0("Baseline gamma: ", x),
      Tlen = function(x) paste0("Time series length: ", x)
    )
  ) +
  theme_bw()



