### Check af remainders ### 

## vi skal plotte hvordan de går mod 0. 
## vi vil gerne overveje hvordan cross-fitting spiller en rolle

source("estimate_test.R")
source("simulate_AR_process.R")
source("conditional_distributions.R")

#true conditional characteristic functions
remainder_true_ccfs = list(
  true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
  true_psi = function(x, u, A, t) char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
)

# parameter setup
Ls <- c(1, 10) #without or with cross-fit
B <- 10 #number of evaluation points
Tlens <- c(20, 50, 100, 200, 500, 1000, 1500, 2000, 5000, 10000) #different chain lengths
Tmax <- max(Tlens) #maximum chain length
M <- 200 #Monte Carlo samples

#generate matrix for VAR(1)
set.seed(420)
A <- matrix(runif(16, -1, 1), nrow = 4) %>% round(2)
A[1, -1] <- 0
d <- ncol(A) - 1

# simulate mu and nu
mu <- matrix(rnorm(B), ncol = 1)
nu <- matrix(rnorm(B * d), ncol = d)


#make grid for results
n_rows <- M * length(Tlens) * length(Ls)
out    <- vector("list", n_rows)
row    <- 1L

for (m in 1:M) {
  # simulate one long path per replication
  full_data <- simulate_AR_process(Tmax, A = A, d = d)
  
  for (T in Tlens) {
    data_sub <- full_data[1:T, ]
    
    for (L in Ls) {
      
      est <- estimate_stat(
        data = data_sub, 
        L = L, 
        B = B,
        mu = mu, nu = nu, 
        A = A, 
        parametric_plugin_AR1 = FALSE,
        remainder_true_ccfs = remainder_true_ccfs
      )
      
      out[[row]] <- tibble(
        Tlen = T,
        L    = L,
        rep  = m,
        R1   = est$Remainders$R1,
        R2   = est$Remainders$R2,
        R3   = est$Remainders$R3,
        S_hat = est$S_hat,
        S = est$S_true
      )
      row <- row + 1L
    }
  }
  
  message("repetition = ", m)
}

results <- bind_rows(out)

saveRDS(results, 
        file = "datasets/remainders/remainders_w_wo_crossfit_B_10_Tlen_20_50_100_200_500_1000_1500_2000_5000_10000.rds")

remainder_df <- readRDS("datasets/remainders/remainders_w_wo_crossfit_B_10_Tlen_20_50_100_200_500_1000_1500_2000_5000_10000.rds")

remainder_df_long <- remainder_df %>%
  pivot_longer(
    cols = c(R1, R2, R3),
    names_to = "remainder",
    values_to = "value"
  ) %>% 
  filter(L == 10) #only look at cross-fitting

remainder_df_long <- remainder_df_long %>% 
  mutate(scaled_value = sqrt(log(Tlen)) * value,
         Tlen_f = factor(Tlen))

initial_values <- remainder_df_long %>% 
  filter(Tlen == 20) %>% 
  dplyr::select(-c(Tlen, S, S_hat)) %>% 
  rename(value_init = value) 
  
remainder_same_scale_df <- remainder_df_long %>% 
  full_join(initial_values, by = c("rep", "L", "remainder")) %>% 
  mutate(value = value/value_init)



remainder_summary <- remainder_df_long %>%
  group_by(Tlen, L, remainder) %>%
  summarize(
    mean_value = mean(scaled_value),
    .groups = "drop") %>% 
  mutate(Tlen_f = factor(Tlen))

ggplot() +
  # violins of blown-up remainders at each Tlen
  geom_violin(
    data = remainder_df_long,
    aes(x = Tlen, y = scaled_value, group = Tlen_f),
    fill = "grey80",
    color = "grey60",
    alpha = 0.7
  ) +
  # mean curve over Tlen (inside each facet)
  geom_line(
    data = remainder_summary,
    aes(x = Tlen, y = mean_value, group = 1),
    color = "red",
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = remainder_summary,
    aes(x = Tlen, y = mean_value),
    color = "red",
    size = 1.8,
    inherit.aes = FALSE
  ) +
  facet_wrap(~remainder, labeller = label_both,
             scales = "free") +
  labs(
    x = "T (chain length)",
    y = "Scaled remainder: R_k * sqrt(log T)",
    title = "Blown-up remainders with violin distributions over T",
    subtitle = "Violin = distribution over reps; red line = Monte Carlo mean"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


ggplot() +
  # all individual rep curves (grey, transparent)
  geom_line(
    data = remainder_df_long,
    aes(x = Tlen, y = scaled_value, group = interaction(rep, remainder)),
    alpha = 0.5,
    color = "grey50"
  ) +
  # mean over reps, colored by remainder
  geom_line(
    data = remainder_summary,
    aes(x = Tlen, y = mean_value, color = remainder),
    linewidth = 0.5
  ) +
  scale_x_log10() +  # optional but nice here
  facet_wrap(~ L, labeller = label_both) +
  labs(
    x = "T (chain length)",
    y = "Remainder value",
    color = "Remainder",
    title = "Remainders R1, R2, R3 vs chain length",
    # subtitle = "Thin lines: individual reps, thick lines: Monte Carlo mean"
  ) +
  theme_bw()




remainder_df_long %>% filter(Tlen == 50, L == 1) %>% 
  ggplot(aes(x = scaled_value)) + 
  # xlim(0, 1) + 
  geom_histogram(color = "white") + 
  facet_wrap(~remainder + Tlen, 
             labeller = labeller(Tlen = function(x) paste0("Tlen = ", x)))
