### Plot document for VAR ###

## Cross fitting ## --- perhaps increasing B???
no_cross_fit_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/L_1_Tlen_30_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_1_Tlen_100_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_1_Tlen_500_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_1_Tlen_1000_B_10.rds")
)
cross_fit_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/L_10_Tlen_30_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_100_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_500_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_B_10.rds")
)

assess_cross_fit_df <- dplyr::bind_rows(
  no_cross_fit_df %>% dplyr::mutate(crossfit = "Without cross-fitting"),
  cross_fit_df %>% dplyr::mutate(crossfit = "With cross-fitting")
)

asses_cross_fit_plot <- assess_cross_fit_df %>% 
  filter(Tlen != 30) %>% 
  dplyr::select(-starts_with("se"), -rate_parametric) %>% 
  tidyr::pivot_longer(
    cols = c(rate_nonparametric),
    values_to = "rate",
    names_to  = "method",
    names_prefix = "rate_"
  ) %>% 
  ggplot(aes(x = alpha, y = rate)) +
  geom_line() +
  geom_abline(color = "darkgrey", linetype = "dashed") +
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate"
  ) +
  facet_grid(
    crossfit ~ Tlen + B + baseline_gamma,
    labeller = label_bquote(
      cols = T == .(Tlen) ~ "," ~ B == .(B)
    )
  ) +
  coord_fixed() +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    legend.position  = "bottom"
  )

ggsave("images/cross-fitting.eps",
       width = 500,
       height = 370,
       units = "px",
       scale = 4,
       plot = asses_cross_fit_plot)



## fixed B ##
level_fixed_B_VAR_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/L_10_Tlen_100_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_500_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_B_10.rds")
)

level_fixed_B_plot <-  level_fixed_B_VAR_df %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
               values_to = "rate",
               names_to = "Method",
               names_prefix = "rate_") %>% 
  mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
  mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
  mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
  ggplot(aes(x = alpha, y = rate, color = Method)) + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate"
  ) + 
  geom_line() + 
  coord_fixed() + 
  geom_abline(color = "grey", linetype = "dashed") + 
  facet_wrap(~Tlen,
             labeller = label_bquote(
               cols = T == .(Tlen)
             )) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom")

ggsave("images/level_fixed_B_VAR.eps",
       width = 500,
       height = 250,
       units = "px",
       scale = 4,
       plot = level_fixed_B_plot)




# Local alternatives #

local_alternative_VAR_fixed_B_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/_Tlen_100_baseline_gamma_3_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_100_baseline_gamma_5_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_3_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_5_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_10_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_15_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_3_B_10.rds"),
  readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_5_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_10_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_15_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_3_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_5_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_10_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_15_B_10.rds")
)

#full plot for local alternative
local_alternative_VAR_fixed_B_full_plot <- local_alternative_VAR_fixed_B_df %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
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
  facet_grid(baseline_gamma ~ Tlen,
             labeller = label_bquote(
               rows = gamma[0] == .(baseline_gamma), 
               cols = Tlen   == .(Tlen)             
             )) + 
  theme_bw() + 
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    ),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom")
local_alternative_VAR_fixed_B_full_plot

## height, width find fra kalibrering med store T ##
ggsave("images/VAR_power_full_plot.eps",
       width = 500,
       height = 450,
       units = "px",
       scale = 4,
       plot = local_alternative_VAR_fixed_B_full_plot)

local_alternative_VAR_fixed_B_alpha_0.05_plot <- local_alternative_VAR_fixed_B_df %>% 
  filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
               values_to = "rate",
               names_to = "Method",
               names_prefix = "rate_") %>% 
  mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
  mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
  mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
  ggplot(aes(x = Tlen, y = rate, color = Method)) + 
  geom_line() + 
  labs(x = "Tlen",
       y = expression(paste("Estimated rejection rate (", alpha, " = 5%)"))) + 
  ylim(0, 1) +
  geom_point(size = 0.5) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~baseline_gamma + B,
             labeller = label_bquote(
               cols = gamma[0] == .(baseline_gamma) ~ "," ~ B == .(B)
             )) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom")
local_alternative_VAR_fixed_B_alpha_0.05_plot


## Increasing B ##

# Level #
level_increasing_B_VAR_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/_Tlen_100_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_0_B_22.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_0_B_30.rds")
  # readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_0_B_42.rds")
)

level_increasing_B_plot <-  level_increasing_B_VAR_df %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
               values_to = "rate",
               names_to = "Method",
               names_prefix = "rate_") %>% 
  mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
  mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
  mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
  ggplot(aes(x = alpha, y = rate, color = Method)) + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate"
  ) + 
  geom_line() + 
  coord_fixed() + 
  geom_abline(color = "grey", linetype = "dashed") + 
  facet_wrap(~Tlen + B,
             labeller = label_bquote(
               cols = T == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NA),
        legend.position = "bottom")

ggsave("images/level_increasing_B_VAR.eps",
       width = 500,
       height = 250,
       units = "px",
       scale = 4,
       plot = level_increasing_B_plot)


# Local alternatives #
local_alternative_VAR_increasing_B_df <- comb_rej_rate_large_obj(
  readRDS("datasets/new_sims/VAR/_Tlen_100_baseline_gamma_3_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_100_baseline_gamma_5_B_10.rds"),
  # readRDS("datasets/new_sims/VAR/_Tlen_100_baseline_gamma_10_B_10.rds"),
  # readRDS("datasets/new_sims/VAR/_Tlen_100_baseline_gamma_15_B_10.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_3_B_22.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_5_B_22.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_10_B_22.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_15_B_22.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_3_B_30.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_5_B_30.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_10_B_30.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_15_B_30.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_3_B_42.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_5_B_42.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_10_B_42.rds"),
  readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_15_B_42.rds")
# 
#   readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_3_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/L_10_Tlen_500_baseline_gamma_5_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_10_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_500_baseline_gamma_15_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_3_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/L_10_Tlen_1000_baseline_gamma_5_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_10_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_1000_baseline_gamma_15_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_3_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_5_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_10_B_10.rds"),
#   readRDS("datasets/new_sims/VAR/_Tlen_2000_baseline_gamma_15_B_10.rds")
)

#full plot for local alternative
local_alternative_VAR_increasing_B_full_plot <- local_alternative_VAR_increasing_B_df %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
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
local_alternative_VAR_increasing_B_full_plot

local_alternative_VAR_increasing_B_alpha_0.05_plot <- local_alternative_VAR_increasing_B_df %>% 
  filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
               values_to = "rate",
               names_to = "Method",
               names_prefix = "rate_") %>% 
  mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
  mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
  mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
  ggplot(aes(x = Tlen, y = rate, color = Method)) + 
  geom_line() + 
  labs(x = "T",
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
local_alternative_VAR_increasing_B_alpha_0.05_plot

ggsave("images/local_alternatives_increasing_B_VAR.eps",
       width = 1000,
       height = 500,
       units = "px",
       scale = 2.5,
       plot = local_alternative_VAR_increasing_B_alpha_0.05_plot)

## combining fixed B and increasing B local alternatives ##
local_alt_all_B <- bind_rows(
  local_alternative_VAR_increasing_B_df %>% 
    filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
    pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
                 values_to = "rate",
                 names_to = "Method",
                 names_prefix = "rate_") %>% 
    mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
    mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
    mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
    mutate(B_scheme = "increasing"),
  local_alternative_VAR_fixed_B_df %>% 
    filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
    pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
                 values_to = "rate",
                 names_to = "Method",
                 names_prefix = "rate_") %>% 
    mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
    mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
    mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
    mutate(B_scheme = "fixed")
) %>% mutate(B_scheme = factor(B_scheme, levels = c("increasing", "fixed")))

local_alternatives_combined_plot <- local_alt_all_B %>% 
  ggplot(aes(x = Tlen, y = rate, color = Method, linetype = B_scheme)) + 
  geom_line() + 
  labs(x = "T",
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
        legend.position = "bottom",
        legend.box = "vertical",
        legend.spacing.y = unit(-0.5, "lines"),
        legend.box.spacing = unit(0, "lines")) + 
  scale_linetype_manual(
    name = expression(B(T)),
    values = c("increasing" = "solid", "fixed" = "dashed"),
    labels = c(latex2exp::TeX('$\\lceil T^{49/100} \\rceil$'), 10)
  )
local_alternatives_combined_plot


ggsave("images/local_alternatives_combined_B.eps",
       width = 500,
       height = 350,
       units = "px",
       scale = 4,
       plot = local_alternatives_combined_plot)



local_alt_all_B2 <- bind_rows(
  local_alternative_VAR_increasing_B_df %>% 
    # filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
    pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
                 values_to = "rate",
                 names_to = "Method",
                 names_prefix = "rate_") %>% 
    mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
    mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
    mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
    mutate(B_scheme = "increasing"),
  local_alternative_VAR_fixed_B_df %>% 
    # filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
    pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin),
                 values_to = "rate",
                 names_to = "Method",
                 names_prefix = "rate_") %>% 
    mutate(Method = ifelse(Method == "oracle_plugin", "Oracle", Method)) %>% 
    mutate(Method = ifelse(Method == "parametric_plugin", "Parametric", Method)) %>% 
    mutate(Method = ifelse(Method == "nonparametric", "Nonparametric", Method)) %>% 
    mutate(B_scheme = "fixed")
) %>% mutate(B_scheme = factor(B_scheme, levels = c("increasing", "fixed")))

local_alternatives_combined_plot2 <- local_alt_all_B2 %>% 
  ggplot(aes(x = alpha, y = rate, color = Method, linetype = B_scheme)) + 
  geom_line() + 
  coord_fixed() + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Rejection rate") +
  ylim(0, 1) + 
  facet_grid(baseline_gamma ~ Tlen,
             labeller = label_bquote(
               rows = gamma[0] == .(baseline_gamma), 
               cols = T == .(Tlen)  
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(-0.5, "lines"),
    legend.box.spacing = unit(0, "lines"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      vjust = 1
    )) + 
  scale_linetype_manual(
    name = expression(B(T)),
    values = c("increasing" = "solid", "fixed" = "dashed"),
    labels = c(latex2exp::TeX('$\\lceil T^{49/100} \\rceil$'), 10)
  )

ggsave("images/VAR_full_plot_local_alternatives.eps",
       width = 500,
       height = 550,
       units = "px",
       scale = 4,
       local_alternatives_combined_plot2)
  
  # ggplot(aes(x = alpha, y = rate, color = Method, linetype = B_scheme)) + 
  # geom_line() + 
  # labs(x = "T",
  #      y = expression(paste("Estimated rejection rate (", alpha, " = 5%)"))) + 
  # ylim(0, 1) +
  # geom_point(size = 0.5) + 
  # geom_abline(color = "darkgrey", linetype = "dashed") + 
  # facet_wrap(~baseline_gamma,
  #            labeller = label_bquote(
  #              cols = gamma[0] == .(baseline_gamma) 
  #            )) + 
  # theme_bw() + 
  # theme(strip.background = element_rect(fill = "white", color = NA),
  #       legend.position = "bottom",
  #       legend.box = "vertical",
  #       legend.spacing.y = unit(-0.5, "lines"),
  #       legend.box.spacing = unit(0, "lines")) + 
  # scale_linetype_manual(
  #   name = expression(B(T)),
  #   values = c("increasing" = "solid", "fixed" = "dashed"),
  #   labels = c(latex2exp::TeX('$\\lceil T^{49/100} \\rceil$'), 10)
  # )

