library(patchwork)
library(tidyverse)
rm(list = ls())

#T1000
T1000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_1.rds")
T1000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_5.rds")
T1000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_9.rds")

#T2000
T2000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_1.rds")
T2000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_5.rds")
T2000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_9.rds")

#T3000
T3000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_1.rds")
T3000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_5.rds")
T3000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_9.rds")

#T4000
T4000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_1.rds")
T4000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_5.rds")
T4000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_9.rds")

#T5000
T5000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_1.rds")
T5000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_5.rds")
T5000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_9.rds")


T1000_bgamm1$sim_rej_obj$rejection_rate_df %>% 
  ggplot(aes(x = alpha, y = rate_parametric)) + 
  geom_line() + 
  geom_abline(color = "red")


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

test <- comb_rej_rate_large_obj(T1000_bgamm1, T1000_bgamm5, T1000_bgamm9, 
                                T2000_bgamm1, T2000_bgamm5, T2000_bgamm9, 
                                T3000_bgamm1, T3000_bgamm5, T3000_bgamm9, 
                                T4000_bgamm1, T4000_bgamm5, T4000_bgamm9, 
                                T5000_bgamm1, T5000_bgamm5, T5000_bgamm9)

#make dataframe long for pivoting
test_long <- test %>% 
  dplyr::select(-starts_with("se_")) %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin), 
               values_to = "rate", 
               names_to = "method",
               names_prefix = "rate_")


ggplot(test_long, aes(x = alpha, y = rate, color = method)) +
  geom_line() +
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") + 
  geom_abline(color = "darkgrey", linetype = "dashed") +
  facet_grid(
    baseline_gamma ~ Tlen,
    labeller = label_bquote(
      rows = gamma[0] == .(baseline_gamma),  # row labels: γ₀ = value
      cols = T[len]   == .(Tlen)             # col labels: T_len = value
    )
  ) +
  theme_bw() + 
  theme(
    strip.background = element_rect(
      fill  = "white",   # or "grey90", "#f0f0f0", etc.
      colour = NA        # or "black" if you want a border
    ),
    legend.position = "bottom"
  )
  


plot_rej_rate_large_obj <- function(large_obj) {
  metadata <- large_obj$metadata
  
  title <- paste0("Chain length: ", metadata$Tlen, ", ", 
                  "baseline gamma: ", metadata$baseline_gamma, ", ",  
                  "actual gamma: ", round(metadata$actual_gamma, 2))
  
  subtitle <- paste0("B: ", metadata$B, ", ",
                     "L: ", metadata$L, ", ",
                     "repetitions: ", metadata$repetitions)
  
  rej_rate_df <- large_obj$sim_rej_obj$rejection_rate_df
  
  rej_rate_df %>% 
    ggplot(aes(x = alpha, y = rate_nonparametric)) + 
    geom_line() + 
    geom_abline(color = "red") +
    ylim(0, 1) + 
    xlab(expression(paste("Significance level (", alpha, ")"))) + 
    ylab("Est. rejection rate") + 
    ggtitle(title, subtitle = subtitle) +
    theme_bw()
}

#T1000
p1000_1 <- plot_rej_rate_large_obj(T1000_bgamm1)
p1000_5 <- plot_rej_rate_large_obj(T1000_bgamm5)

#T2000
p2000_1 <- plot_rej_rate_large_obj(T2000_bgamm1)
p2000_5 <- plot_rej_rate_large_obj(T2000_bgamm5)
plot_rej_rate_large_obj(T2000_bgamm9)

#T3000
plot_rej_rate_large_obj(T3000_bgamm1)

#T5000
plot_rej_rate_large_obj(T5000_bgamm1)


#T5000 calibration
T5000_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_0.rds")
plot_rej_rate_large_obj(T5000_bgamm0)


### --------- calibration plots --------- ###

# AR1 #
T100_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_100_baseline_gamma_0.rds")
T200_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_200_baseline_gamma_0.rds")
T300_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_300_baseline_gamma_0.rds")
T400_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_400_baseline_gamma_0.rds")
T500_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_500_baseline_gamma_0.rds")
T600_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_600_baseline_gamma_0.rds")
T700_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_700_baseline_gamma_0.rds")
T800_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_800_baseline_gamma_0.rds")
T900_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_900_baseline_gamma_0.rds")
T1000_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_0.rds")
T2000_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_0.rds")
T3000_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_0.rds")
T4000_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_0.rds")
T5000_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_0.rds")
T10000_bgamm0 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_10000_baseline_gamma_0.rds")

calibration_big_df <- comb_rej_rate_large_obj(
  T100_bgamm0,
  T200_bgamm0,
  T300_bgamm0,
  T400_bgamm0,
  T500_bgamm0,
  T600_bgamm0,
  T700_bgamm0,
  T800_bgamm0,
  T900_bgamm0,
  T1000_bgamm0,
  T2000_bgamm0,
  T3000_bgamm0,
  T4000_bgamm0,
  T5000_bgamm0,
  T10000_bgamm0
)

calibration_big_df_long <- calibration_big_df %>% 
  dplyr::select(-starts_with("se_")) %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin), 
               values_to = "rate", 
               names_to = "method",
               names_prefix = "rate_")

ggplot(calibration_big_df_long, aes(x = alpha, y = rate, color = method)) +
  geom_line() +
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") + 
  geom_abline(color = "darkgrey", linetype = "dashed") +
  facet_wrap(~Tlen + B + baseline_gamma,
             ncol = 5,
             labeller = label_bquote(
               Tlen == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  # facet_grid( #maybe make face wrap because only one gamma
  #   rows = vars(baseline_gamma),
  #   cols = vars(Tlen, B),
  #   # baseline_gamma ~ Tlen + B,
  #   labeller = label_bquote(
  #     rows = gamma[0] == .(baseline_gamma),  # row labels: γ₀ = value
  #     cols = Tlen ==  .(Tlen) ~ ","~ B == .(B)  # col labels: T_len = value
  #   ),
  # ) +
  theme_bw() + 
  theme(
    strip.background = element_rect(
      fill  = "white",   # or "grey90", "#f0f0f0", etc.
      colour = NA        # or "black" if you want a border
    ),
    legend.position = "bottom"
  )

### non parametric ###
calibration_big_df_long_nonpar <- calibration_big_df_long %>% 
  dplyr::filter(method == "nonparametric")

ggplot(calibration_big_df_long_nonpar, aes(x = alpha, y = rate)) +
  geom_line() +
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") + 
  geom_abline(color = "darkgrey", linetype = "dashed") +
  facet_wrap(~Tlen + B + baseline_gamma,
             ncol = 5,
             labeller = label_bquote(
               Tlen == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  # facet_grid( #maybe make face wrap because only one gamma
  #   rows = vars(baseline_gamma),
  #   cols = vars(Tlen, B),
  #   # baseline_gamma ~ Tlen + B,
  #   labeller = label_bquote(
  #     rows = gamma[0] == .(baseline_gamma),  # row labels: γ₀ = value
  #     cols = Tlen ==  .(Tlen) ~ ","~ B == .(B)  # col labels: T_len = value
  #   ),
  # ) +
  theme_bw() + 
  theme(
    strip.background = element_rect(
      fill  = "white",   # or "grey90", "#f0f0f0", etc.
      colour = NA        # or "black" if you want a border
    ),
    legend.position = "bottom"
  )


# alpha 5% table

tab <- calibration_big_df_long %>%
  filter(0.0475 <= alpha,
         alpha <= 0.0525,
         method == "nonparametric") %>%
  arrange(Tlen) %>%
  dplyr::select(Tlen, rate)

n  <- nrow(tab)
n1 <- ceiling(n / 2)

tab1 <- tab[1:n1, ]
tab2 <- tab[(n1 + 1):n, ]

tab_wide <- dplyr::tibble(
  Tlen_1 = tab1$Tlen,
  rate_1 = tab1$rate,
  Tlen_2 = c(tab2$Tlen, rep(NA, n1 - nrow(tab2))),
  rate_2 = c(tab2$rate, rep(NA, n1 - nrow(tab2)))
)

tab_wide %>%
  kableExtra::kbl(
    format    = "latex",
    booktabs  = TRUE,
    col.names = c("$T_{\\text{len}}$", "$\\hat r$", "$T_{\\text{len}}$", "$\\hat r$"),
    digits    = 3,
    escape    = FALSE,
    na        = ""
  ) %>%
  kableExtra::kable_styling(
    position   = "center",
    full_width = FALSE
  )

# CIR #
T100_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_100_baseline_gamma_0.rds")
T200_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_200_baseline_gamma_0.rds")
T300_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_300_baseline_gamma_0.rds")
T400_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_400_baseline_gamma_0.rds")
T500_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_500_baseline_gamma_0.rds")
T600_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_600_baseline_gamma_0.rds")
T700_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_700_baseline_gamma_0.rds")
T800_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_800_baseline_gamma_0.rds")
T900_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_900_baseline_gamma_0.rds")
T1000_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_1000_baseline_gamma_0.rds")
T2000_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_2000_baseline_gamma_0.rds")
T3000_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_3000_baseline_gamma_0.rds")
T4000_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_4000_baseline_gamma_0.rds")
T5000_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_5000_baseline_gamma_0.rds")
T10000_bgamm0_CIR <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_CIR_4D_local_alt_Tlen_10000_baseline_gamma_0.rds")

calibration_big_df_CIR <- comb_rej_rate_large_obj(
  T100_bgamm0_CIR,
  T200_bgamm0_CIR,
  T300_bgamm0_CIR,
  T400_bgamm0_CIR,
  T500_bgamm0_CIR,
  T600_bgamm0_CIR,
  T700_bgamm0_CIR,
  T800_bgamm0_CIR,
  T900_bgamm0_CIR,
  T1000_bgamm0_CIR,
  T2000_bgamm0_CIR,
  T3000_bgamm0_CIR,
  T4000_bgamm0_CIR,
  T5000_bgamm0_CIR,
  T10000_bgamm0_CIR
)

calibration_big_df_CIR_long <- calibration_big_df_CIR %>% 
  dplyr::select(-starts_with("se_")) %>% 
  pivot_longer(cols = c(rate_nonparametric), 
               values_to = "rate", 
               names_to = "method",
               names_prefix = "rate_")

ggplot(calibration_big_df_CIR_long, aes(x = alpha, y = rate)) +
  geom_line() +
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") + 
  geom_abline(color = "darkgrey", linetype = "dashed") +
  facet_wrap(~Tlen + B + baseline_gamma,
             ncol = 5,
             labeller = label_bquote(
               Tlen == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  # facet_grid( #maybe make face wrap because only one gamma
  #   rows = vars(baseline_gamma),
  #   cols = vars(Tlen, B),
  #   # baseline_gamma ~ Tlen + B,
  #   labeller = label_bquote(
  #     rows = gamma[0] == .(baseline_gamma),  # row labels: γ₀ = value
  #     cols = Tlen ==  .(Tlen) ~ ","~ B == .(B)  # col labels: T_len = value
  #   ),
  # ) +
  theme_bw() + 
  theme(
    strip.background = element_rect(
      fill  = "white",   # or "grey90", "#f0f0f0", etc.
      colour = NA        # or "black" if you want a border
    ),
    legend.position = "bottom"
  )

# alpha 5% table

tab_CIR <- calibration_big_df_CIR_long %>%
  filter(0.0475 <= alpha,
         alpha <= 0.0525,
         method == "nonparametric") %>%
  arrange(Tlen) %>%
  dplyr::select(Tlen, rate)

n_CIR  <- nrow(tab_CIR)
n1_CIR <- ceiling(n_CIR / 2)

tab1_CIR <- tab_CIR[1:n1_CIR, ]
tab2_CIR <- tab_CIR[(n1_CIR + 1):n_CIR, ]

tab_wide_CIR <- dplyr::tibble(
  Tlen_1 = tab1_CIR$Tlen,
  rate_1 = tab1_CIR$rate,
  Tlen_2 = c(tab2_CIR$Tlen, rep(NA, n1_CIR - nrow(tab2_CIR))),
  rate_2 = c(tab2_CIR$rate, rep(NA, n1_CIR - nrow(tab2_CIR)))
)

tab_wide_CIR %>%
  kableExtra::kbl(
    format    = "latex",
    booktabs  = TRUE,
    col.names = c("$T_{\\text{len}}$", "$\\hat r$", "$T_{\\text{len}}$", "$\\hat r$"),
    digits    = 3,
    escape    = FALSE,
    na        = ""
  ) %>%
  kableExtra::kable_styling(
    position   = "center",
    full_width = FALSE
  )



### ------ power things ------ ###

# AR1 #
#T1000
T1000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_1.rds")
T1000_bgamm3 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_3.rds")
T1000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_5.rds")
T1000_bgamm7 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_7.rds")
T1000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_9.rds")
T1000_bgamm11 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_11.rds")
T1000_bgamm13 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_13.rds")
T1000_bgamm15 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_15.rds")
T1000_bgamm17 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_1000_baseline_gamma_17.rds")

#T2000
T2000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_1.rds")
T2000_bgamm3 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_3.rds")
T2000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_5.rds")
T2000_bgamm7 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_7.rds")
T2000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_9.rds")
T2000_bgamm11 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_11.rds")
T2000_bgamm13 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_13.rds")
T2000_bgamm15 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_15.rds")
T2000_bgamm17 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_2000_baseline_gamma_17.rds")

#T3000
T3000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_1.rds")
T3000_bgamm3 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_3.rds")
T3000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_5.rds")
T3000_bgamm7 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_7.rds")
T3000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_9.rds")
T3000_bgamm11 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_11.rds")
T3000_bgamm13 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_13.rds")
T3000_bgamm15 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_15.rds")
T3000_bgamm17 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_3000_baseline_gamma_17.rds")

#T4000
T4000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_1.rds")
T4000_bgamm3 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_3.rds")
T4000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_5.rds")
T4000_bgamm7 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_7.rds")
T4000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_9.rds")
T4000_bgamm11 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_11.rds")
T4000_bgamm13 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_13.rds")
T4000_bgamm15 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_15.rds")
T4000_bgamm17 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_4000_baseline_gamma_17.rds")

#T5000
T5000_bgamm1 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_1.rds")
T5000_bgamm3 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_3.rds")
T5000_bgamm5 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_5.rds")
T5000_bgamm7 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_7.rds")
T5000_bgamm9 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_9.rds")
T5000_bgamm11 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_11.rds")
T5000_bgamm13 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_13.rds")
T5000_bgamm15 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_15.rds")
T5000_bgamm17 <- readRDS("datasets/local_alternatives_4D/sim_rej_rate_AR1_4D_local_alt_Tlen_5000_baseline_gamma_17.rds")



power_df <- comb_rej_rate_large_obj(T1000_bgamm1, T1000_bgamm3, T1000_bgamm5, T1000_bgamm7, T1000_bgamm9, T1000_bgamm9, T1000_bgamm11, T1000_bgamm13, T1000_bgamm15,
                                    T2000_bgamm1, T2000_bgamm3, T2000_bgamm5, T2000_bgamm7, T2000_bgamm9, T2000_bgamm9, T2000_bgamm11, T2000_bgamm13, T2000_bgamm15,
                                    T3000_bgamm1, T3000_bgamm3, T3000_bgamm5, T3000_bgamm7, T3000_bgamm9, T3000_bgamm9, T3000_bgamm11, T3000_bgamm13, T3000_bgamm15,
                                    T4000_bgamm1, T4000_bgamm3, T4000_bgamm5, T4000_bgamm7, T4000_bgamm9, T4000_bgamm9, T4000_bgamm11, T4000_bgamm13, T4000_bgamm15,
                                    T5000_bgamm1, T5000_bgamm3, T5000_bgamm5, T5000_bgamm7, T5000_bgamm9, T5000_bgamm9, T5000_bgamm11, T5000_bgamm13, T5000_bgamm15)

#only pick alpha = 0.05 / 5%
power_df_alpha_0.05 <- power_df %>% filter(0.0475 <= alpha & alpha <= 0.0525)

power_df_alpha_0.05_long <- power_df_alpha_0.05 %>% 
  dplyr::select(-starts_with("se_")) %>% 
  pivot_longer(cols = c(rate_nonparametric, rate_parametric_plugin, rate_oracle_plugin), 
               values_to = "rate", 
               names_to = "method",
               names_prefix = "rate_")

power_df_alpha_0.05_long_nonparametric <- power_df_alpha_0.05_long %>% 
  dplyr::filter(method == "nonparametric")

power_df_alpha_0.05_long_nonparametric %>% 
  # ggplot(aes(x = Tlen, y = rate, color = method)) + 
  ggplot(aes(x = Tlen, y = rate)) + 
  geom_line() + 
  labs(x = "Tlen",
       y = expression(paste("Est. rejection rate (", alpha, " = 5%)"))) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~baseline_gamma,
             labeller = label_bquote(
               rows = gamma[0] == .(baseline_gamma)  # row labels: γ₀ = value 
             )) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.background = element_rect(
          fill = "white", color = NA
        ))


# CIR #
T1000_bgamm1_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_1000_baseline_gamma_1.rds")
T1000_bgamm3_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_1000_baseline_gamma_3.rds")
T1000_bgamm5_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_1000_baseline_gamma_5.rds")
T1000_bgamm7_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_1000_baseline_gamma_7.rds")
T1000_bgamm9_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_1000_baseline_gamma_9.rds")
T1000_bgamm11_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_1000_baseline_gamma_11.rds")

#T2000
T2000_bgamm1_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_2000_baseline_gamma_1.rds")
T2000_bgamm3_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_2000_baseline_gamma_3.rds")
T2000_bgamm5_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_2000_baseline_gamma_5.rds")
T2000_bgamm7_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_2000_baseline_gamma_7.rds")
T2000_bgamm9_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_2000_baseline_gamma_9.rds")
T2000_bgamm11_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_2000_baseline_gamma_11.rds")

#T3000
T3000_bgamm1_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_3000_baseline_gamma_1.rds")
T3000_bgamm3_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_3000_baseline_gamma_3.rds")
T3000_bgamm5_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_3000_baseline_gamma_5.rds")
T3000_bgamm7_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_3000_baseline_gamma_7.rds")
T3000_bgamm9_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_3000_baseline_gamma_9.rds")
T3000_bgamm11_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_3000_baseline_gamma_11.rds")

#T4000
T4000_bgamm1_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_4000_baseline_gamma_1.rds")
T4000_bgamm3_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_4000_baseline_gamma_3.rds")
T4000_bgamm5_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_4000_baseline_gamma_5.rds")
T4000_bgamm7_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_4000_baseline_gamma_7.rds")
T4000_bgamm9_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_4000_baseline_gamma_9.rds")
T4000_bgamm11_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_4000_baseline_gamma_11.rds")

#T5000
T5000_bgamm1_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_5000_baseline_gamma_1.rds")
T5000_bgamm3_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_5000_baseline_gamma_3.rds")
T5000_bgamm5_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_5000_baseline_gamma_5.rds")
T5000_bgamm7_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_5000_baseline_gamma_7.rds")
T5000_bgamm9_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_5000_baseline_gamma_9.rds")
T5000_bgamm11_CIR <- readRDS("datasets/local_alternatives_4D/2nd_sim_rej_rate_CIR_4D_local_alt_Tlen_5000_baseline_gamma_11.rds")

power_df_CIR <- comb_rej_rate_large_obj(T1000_bgamm1_CIR, T1000_bgamm3_CIR, T1000_bgamm5_CIR, T1000_bgamm7_CIR, T1000_bgamm9_CIR, T1000_bgamm9_CIR, T1000_bgamm11_CIR, 
                                    T2000_bgamm1_CIR, T2000_bgamm3_CIR, T2000_bgamm5_CIR, T2000_bgamm7_CIR, T2000_bgamm9_CIR, T2000_bgamm9_CIR, T2000_bgamm11_CIR,
                                    T3000_bgamm1_CIR, T3000_bgamm3_CIR, T3000_bgamm5_CIR, T3000_bgamm7_CIR, T3000_bgamm9_CIR, T3000_bgamm9_CIR, T3000_bgamm11_CIR,
                                    T4000_bgamm1_CIR, T4000_bgamm3_CIR, T4000_bgamm5_CIR, T4000_bgamm7_CIR, T4000_bgamm9_CIR, T4000_bgamm9_CIR, T4000_bgamm11_CIR,
                                    T5000_bgamm1_CIR, T5000_bgamm3_CIR, T5000_bgamm5_CIR, T5000_bgamm7_CIR, T5000_bgamm9_CIR, T5000_bgamm9_CIR, T5000_bgamm11_CIR
)

power_df_CIR_alpha_0.05 <- power_df_CIR %>% filter(0.0475 <= alpha & alpha <= 0.0525)


power_df_CIR_alpha_0.05_long <- power_df_CIR_alpha_0.05 %>% 
  dplyr::select(-starts_with("se_")) %>% 
  pivot_longer(cols = c(rate_nonparametric), 
               values_to = "rate", 
               names_to = "method",
               names_prefix = "rate_")

power_df_CIR_alpha_0.05_long %>% 
  # ggplot(aes(x = Tlen, y = rate, color = method)) + 
  ggplot(aes(x = Tlen, y = rate)) + 
  geom_line() + 
  labs(x = "Tlen",
       y = expression(paste("Est. rejection rate (", alpha, " = 5%)"))) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~baseline_gamma,
             labeller = label_bquote(
               rows = gamma[0] == .(baseline_gamma)  # row labels: γ₀ = value 
             )) + 
  theme_bw() + 
  theme(legend.position = "bottom",
        strip.background = element_rect(
          fill = "white", color = NA
        ))


