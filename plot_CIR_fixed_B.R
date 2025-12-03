### CIR calibration plots for simulations december 2nd to december 3rd ### 

#function to combine several data frames
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
                            actual_gamma = metadata_i$actual_gamma,
                            L = metadata_i$L)
    
    df <- bind_rows(df, df_i)
  }
  
  df
}

#CIR calibration df
CIR_calibration_fixed_B_df <- comb_rej_rate_large_obj(
  readRDS("datasets/niveau/CIR/grid_Tlen_100_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_200_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_300_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_400_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_500_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_600_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_700_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_800_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_900_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_1000_baseline_gamma_0_B_10.rds"),
  readRDS("datasets/niveau/CIR/grid_Tlen_2000_baseline_gamma_0_B_10.rds")
)


#CIR calibration plot
CIR_calibration_fixed_B_plot <- CIR_calibration_fixed_B_df %>% 
  ggplot(aes(x = alpha, y = rate_nonparametric)) +
  geom_line() +
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate"
  ) +
  facet_wrap(~Tlen,
             labeller = label_bquote(
               cols = Tlen == .(Tlen) ~ "," ~ B == .(B)
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", colour = NA),
    legend.position = "bottom")

CIR_calibration_fixed_B_plot

### CIR local alternatives plots for fixed B ###
local_alt_100_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_100_baseline_gamma_1_B_10.rds")
local_alt_100_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_100_baseline_gamma_10_B_10.rds")
local_alt_100_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_100_baseline_gamma_15_B_10.rds")
local_alt_100_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_100_baseline_gamma_20_B_10.rds")
local_alt_100_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_100_baseline_gamma_30_B_10.rds")

local_alt_200_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_200_baseline_gamma_1_B_10.rds")
local_alt_200_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_200_baseline_gamma_10_B_10.rds")
local_alt_200_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_200_baseline_gamma_15_B_10.rds")
local_alt_200_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_200_baseline_gamma_20_B_10.rds")
local_alt_200_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_200_baseline_gamma_30_B_10.rds")

local_alt_300_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_300_baseline_gamma_1_B_10.rds")
local_alt_300_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_300_baseline_gamma_10_B_10.rds")
local_alt_300_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_300_baseline_gamma_15_B_10.rds")
local_alt_300_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_300_baseline_gamma_20_B_10.rds")
local_alt_300_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_300_baseline_gamma_30_B_10.rds")

local_alt_400_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_400_baseline_gamma_1_B_10.rds")
local_alt_400_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_400_baseline_gamma_10_B_10.rds")
local_alt_400_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_400_baseline_gamma_15_B_10.rds")
local_alt_400_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_400_baseline_gamma_20_B_10.rds")
local_alt_400_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_400_baseline_gamma_30_B_10.rds")

local_alt_500_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_500_baseline_gamma_1_B_10.rds")
local_alt_500_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_500_baseline_gamma_10_B_10.rds")
local_alt_500_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_500_baseline_gamma_15_B_10.rds")
local_alt_500_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_500_baseline_gamma_20_B_10.rds")
local_alt_500_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_500_baseline_gamma_30_B_10.rds")

local_alt_600_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_600_baseline_gamma_1_B_10.rds")
local_alt_600_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_600_baseline_gamma_10_B_10.rds")
local_alt_600_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_600_baseline_gamma_15_B_10.rds")
local_alt_600_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_600_baseline_gamma_20_B_10.rds")
local_alt_600_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_600_baseline_gamma_30_B_10.rds")

local_alt_700_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_700_baseline_gamma_1_B_10.rds")
local_alt_700_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_700_baseline_gamma_10_B_10.rds")
local_alt_700_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_700_baseline_gamma_15_B_10.rds")
local_alt_700_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_700_baseline_gamma_20_B_10.rds")
local_alt_700_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_700_baseline_gamma_30_B_10.rds")

local_alt_800_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_800_baseline_gamma_1_B_10.rds")
local_alt_800_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_800_baseline_gamma_10_B_10.rds")
local_alt_800_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_800_baseline_gamma_15_B_10.rds")
local_alt_800_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_800_baseline_gamma_20_B_10.rds")
local_alt_800_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_800_baseline_gamma_30_B_10.rds")

local_alt_2000_gamma_1 <- readRDS("datasets/niveau/CIR/grid_Tlen_2000_baseline_gamma_1_B_10.rds")
local_alt_2000_gamma_10 <- readRDS("datasets/niveau/CIR/grid_Tlen_2000_baseline_gamma_10_B_10.rds")
local_alt_2000_gamma_15 <- readRDS("datasets/niveau/CIR/grid_Tlen_2000_baseline_gamma_15_B_10.rds")
local_alt_2000_gamma_20 <- readRDS("datasets/niveau/CIR/grid_Tlen_2000_baseline_gamma_20_B_10.rds")
local_alt_2000_gamma_30 <- readRDS("datasets/niveau/CIR/grid_Tlen_2000_baseline_gamma_30_B_10.rds")

CIR_local_alternatives_fixed_B_df <- comb_rej_rate_large_obj(
  local_alt_100_gamma_1, 
  local_alt_100_gamma_10,
  local_alt_100_gamma_15,
  local_alt_100_gamma_20,
  local_alt_100_gamma_30, # end of 100
  local_alt_200_gamma_1, 
  local_alt_200_gamma_10,
  local_alt_200_gamma_15,
  local_alt_200_gamma_20,
  local_alt_200_gamma_30, # end of 22
  local_alt_300_gamma_1, 
  local_alt_300_gamma_10,
  local_alt_300_gamma_15,
  local_alt_300_gamma_20,
  local_alt_300_gamma_30, #end of 33
  local_alt_400_gamma_1, 
  local_alt_400_gamma_10,
  local_alt_400_gamma_15,
  local_alt_400_gamma_20,
  local_alt_400_gamma_30, # end of 400
  local_alt_500_gamma_1, 
  local_alt_500_gamma_10,
  local_alt_500_gamma_15,
  local_alt_500_gamma_20,
  local_alt_500_gamma_30, #end of 500
  local_alt_600_gamma_1, 
  local_alt_600_gamma_10,
  local_alt_600_gamma_15,
  local_alt_600_gamma_20,
  local_alt_600_gamma_30, #end og 600
  local_alt_700_gamma_1, 
  local_alt_700_gamma_10,
  local_alt_700_gamma_15,
  local_alt_700_gamma_20,
  local_alt_700_gamma_30, #end of 700
  local_alt_800_gamma_1, 
  local_alt_800_gamma_10,
  local_alt_800_gamma_15,
  local_alt_800_gamma_20,
  local_alt_800_gamma_30, #end of 800
  local_alt_2000_gamma_1,
  local_alt_2000_gamma_10,
  local_alt_2000_gamma_15,
  local_alt_2000_gamma_20,
  local_alt_2000_gamma_30 #end of 2000
)

CIR_local_alter_fixed_B_full_plot <- CIR_local_alternatives_fixed_B_df %>% 
  pivot_longer(cols = c(rate_nonparametric), 
               values_to = "rate",
               names_to = "method") %>% 
  ggplot(aes(x = alpha, y = rate)) + 
  geom_line() + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  labs(
    x = expression(paste("Significance level (", alpha, ")")),
    y = "Est. rejection rate") +
  ylim(0, 1) + 
  facet_grid(baseline_gamma ~ Tlen,
             labeller = label_bquote(
               rows = gamma[0] == .(baseline_gamma), 
               cols = Tlen   == .(Tlen)             
             )) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = "white", color = NA)
  )

CIR_local_alter_fixed_B_full_plot

CIR_local_alter_fixed_B_alpha_0.05_plot <- CIR_local_alternatives_fixed_B_df %>% 
  filter(0.0475 <= alpha & alpha <= 0.0525) %>% 
  ggplot(aes(x = Tlen, y = rate_nonparametric)) + 
  geom_line() + 
  labs(x = "Tlen",
       y = expression(paste("Estimated rejection rate (", alpha, " = 5%)"))) + 
  ylim(0, 1) +
  geom_point(size = 0.5) + 
  geom_abline(color = "darkgrey", linetype = "dashed") + 
  facet_wrap(~baseline_gamma,
             labeller = label_bquote(
               cols = gamma[0] == .(baseline_gamma) ~ "," ~ B == .(B)
             )) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill = "white", color = NA))

CIR_local_alter_fixed_B_alpha_0.05_plot


