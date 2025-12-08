#### Level anayslis for B scheme ####

level_analysis_B_scheme <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/L_10_Tlen_100_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_200_baseline_gamma_0_B_14.rds"),
  #readRDS("datasets/new_sims/VAR/L_10_Tlen_300_baseline_gamma_0_B_17.rds"),
  #readRDS("datasets/new_sims/VAR/L_10_Tlen_400_baseline_gamma_0_B_19.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_0_B_22.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_700_baseline_gamma_0_B_25.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_0_B_30.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_1500_baseline_gamma_0_B_36.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_2000_baseline_gamma_0_B_42.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_3000_baseline_gamma_0_B_51.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_5000_baseline_gamma_0_B_65.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_7500_baseline_gamma_0_B_80.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_10000_baseline_gamma_0_B_92.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_15000_baseline_gamma_0_B_112.rds")
)


level_analysis_B_scheme %>% 
  pivot_longer(
    starts_with("rate_"),
    names_to = "Method",
    values_to = "rate",
    names_prefix = "rate_"
  ) %>% 
  filter(Method != "parametric") %>% 
  filter(Method == "nonparametric") %>%
  ggplot(aes(x = alpha, y = rate)) +
  geom_line() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  facet_wrap(vars(Tlen, B), ncol = 4,
             labeller = label_bquote(
               cols  = T == .(Tlen) ~ "," ~ B == .(B))
             #   labeller(
             #   Tlen = function(x) paste0("T = ", x),
             #   B = function(x) paste0("B = ", x)
             # )
             ) +
  coord_fixed() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white", color  = NA),
    axis.text.x = element_text(
      angle = 30,      # angle in degrees
      hjust = 1,       # horizontal justification
      vjust = 1        # vertical justification
    )
  ) +
  xlab(expression(paste("Significance level (", alpha, ")"))) +
  ylab("Rejection rate")




ggsave("images/Large_level_B_scheme_VAR.eps", 
       width = 500, 
       height = 450,
       units = "px",
       scale = 4
)





