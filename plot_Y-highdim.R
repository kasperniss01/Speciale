### plot Y highdim ###

increasing_dim_Y_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/highdim/_d_8_Tlen_500_B_22_gamma_0.rds"),
  # readRDS("datasets/new_sims/VAR/highdim/_d_8_Tlen_500_B_22_gamma_0.55.rds"),
  readRDS("datasets/new_sims/VAR/highdim/_d_16_Tlen_500_B_22_gamma_0.rds"),
  # readRDS("datasets/new_sims/VAR/highdim/_d_16_Tlen_500_B_22_gamma_0.63.rds"),
  readRDS("datasets/new_sims/VAR/highdim/_d_32_Tlen_500_B_22_gamma_0.rds"),
  # readRDS("datasets/new_sims/VAR/highdim/_d_32_Tlen_500_B_22_gamma_0.71.rds"),
  readRDS("datasets/new_sims/VAR/highdim/_d_64_Tlen_500_B_22_gamma_0.rds"),
  # readRDS("datasets/new_sims/VAR/highdim/_d_64_Tlen_500_B_22_gamma_0.77.rds"),
  readRDS("datasets/new_sims/VAR/highdim/_d_128_Tlen_500_B_22_gamma_0.rds"),
  # readRDS("datasets/new_sims/VAR/highdim/_d_128_Tlen_500_B_22_gamma_0.84.rds")
  readRDS("datasets/new_sims/VAR/highdim/_d_256_Tlen_500_B_22_gamma_0.rds")
) %>% mutate(
  d = rep(c(8, 16, 32, 64, 128, 256), each = 200)
)

#level
increasing_Y_level_plot <- increasing_dim_Y_df %>% 
  filter(actual_gamma == 0) %>% 
  ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate"
  ) + 
  geom_abline(color = "darkgrey", linetype = 2) + 
  facet_wrap(~d,
             labeller = label_bquote(
               cols = d == .(d)
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", color = NA)
  )

ggsave("images/increasing_dim_Y_level_plot.eps",
       width = 500,
       heigh = 350,
       units = "px",
       scale = 4,
       increasing_Y_level_plot)

#alternatives
increasing_dim_Y_df %>% filter(actual_gamma != 0) %>% 
  ggplot(aes(x = alpha, y = rate_nonparametric)) + 
  geom_line() +
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate"
  ) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~d + actual_gamma,
             labeller = label_bquote(
               cols = d == .(d) ~ "," ~ gamma == .(actual_gamma)
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", color = NA)
  )


# niveau plot Y 10D
