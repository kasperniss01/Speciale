### Plot document for CIR - increasing B ###


# level #
level_increasing_B_CIR_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/CIR/_Tlen_100_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_0_B_22.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_0_B_30.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_0_B_42.rds")
)


level_increasing_B_plot <-  level_increasing_B_CIR_df %>% 
  pivot_longer(cols = c(rate_nonparametric),
               values_to = "rate",
               names_to = "Method",
               names_prefix = "rate_") %>% 
  mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
  mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
  mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
  ggplot(aes(x = alpha, y = rate)) + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate"
  ) + 
  geom_line() + 
  coord_fixed() + 
  geom_abline(color = "grey", linetype = "dashed") + 
  facet_wrap(~Tlen + B,
             nrow = 2, 
             labeller = label_bquote(
               cols = T == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom")

ggsave("images/level_increasing_B_CIR.png", 
       width = 350,
       height = 200,
       units = "px",
       scale = 5,
       plot = level_increasing_B_plot)

# Local alternatives #
local_alternative_CIR_increasing_B_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/CIR/_Tlen_100_baseline_gamma_1_B_10.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_100_baseline_gamma_10_B_10.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_100_baseline_gamma_20_B_10.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_100_baseline_gamma_30_B_10.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_1_B_22.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_10_B_22.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_20_B_22.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_30_B_22.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_1_B_30.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_10_B_30.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_20_B_30.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_30_B_30.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_1_B_42.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_10_B_42.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_20_B_42.rds"),
  readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_30_B_42.rds")
  # 
  #   readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_3_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/L_10_Tlen_500_baseline_gamma_5_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_10_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_500_baseline_gamma_15_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_3_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/L_10_Tlen_1000_baseline_gamma_5_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_10_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_1000_baseline_gamma_15_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_3_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_5_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_10_B_10.rds"),
  #   readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_15_B_10.rds")
)

#full plot for local alternative
local_alternative_CIR_increasing_B_full_plot <- local_alternative_CIR_increasing_B_df %>% 
  pivot_longer(cols = c(rate_nonparametric),
               values_to = "rate",
               names_to = "Method",
               names_prefix = "rate_") %>% 
  mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
  mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
  mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
  ggplot(aes(x = alpha, y = rate, color = Method)) + 
  geom_line() + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate") +
  ylim(0, 1) + 
  facet_grid(baseline_gamma ~ Tlen + B,
             labeller = label_bquote(
               rows = gamma[0] == .(baseline_gamma), 
               cols = T == .(Tlen) ~ "," ~ B == .(B)  
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom")
local_alternative_CIR_increasing_B_full_plot

local_alternative_CIR_increasing_B_alpha_0.05_plot <- local_alternative_CIR_increasing_B_df %>% 
  filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
  pivot_longer(cols = c(rate_nonparametric),
               values_to = "rate",
               names_to = "Method",
               names_prefix = "rate_") %>% 
  mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
  mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
  mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
  ggplot(aes(x = Tlen, y = rate)) + 
  geom_line() + 
  labs(x = "Tlen",
       y = expression(paste("Estimated rejection rate (", alpha, " = 5%)"))) + 
  ylim(0, 1) +
  geom_point(size = 0.5) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~baseline_gamma,
             labeller = label_bquote(
               cols = gamma[0] == .(baseline_gamma) 
             )) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom")
local_alternative_CIR_increasing_B_alpha_0.05_plot

ggsave("images/local_alternatives_CIR.eps",
       width = 500,
       height = 350,
       units = "px",
       scale = 4,
       plot = local_alternative_CIR_increasing_B_alpha_0.05_plot)














readRDS("datasets/new_sims/CIR/_Tlen_2000_baseline_gamma_0_B_42.rds") %>% 
  comb_rej_rate_large_obj() %>% 
  ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() + 
  geom_abline()
