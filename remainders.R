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
  )

remainder_df_long <- remainder_df_long %>% 
  mutate(scaled_value = sqrt(log(Tlen)) * value,
         Tlen_f = factor(Tlen))

remainder_mean <- remainder_df_long %>% 
  group_by(L, remainder, Tlen) %>% 
  summarise(mean_scaled_value = mean(scaled_value),
            mean_value = mean(value))

remainder_labs <- as_labeller(
  c(
    R1 = "max[b%in%group('[',B,']')]~sqrt(T~log(T))~group('|',R[1](b),'|')",
    R2 = "max[b%in%group('[',B,']')]~sqrt(T~log(T))~group('|',R[2](b),'|')",
    R3 = "max[b%in%group('[',B,']')]~sqrt(T~log(T))~group('|',R[3](b),'|')"
  ),
  label_parsed
)

# Row labels:  L = 1, L = 10
L_labs <- function(x) paste0("L = ", x)

#fully scaled remainder
remainder_plot <- ggplot() + 
  # 200 grey lines = individual replications
  geom_line(
    data = remainder_df_long,
    aes(x = Tlen, y = scaled_value, group = rep),
    alpha = 1,
    color = "grey85",
    linewidth = 0.2
  ) +
  # Red line = mean across replications
  geom_line(
    data = remainder_mean,
    aes(x = Tlen, y = mean_scaled_value),
    color = "red",
    linewidth = 0.7
  ) +
  labs(
    x = expression(T),
    y = "Value"
  ) + 
  ylim(0, 1.5) + 
  facet_grid(
    L ~ remainder,
    labeller = labeller(
      L        = L_labs,
      remainder = remainder_labs
    )
  ) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white", color = NA))
remainder_plot

ggsave("images/remainder_plot.eps",
       width = 500,
       height = 300,
       units = "px",
       scale = 4,
       remainder_plot)




