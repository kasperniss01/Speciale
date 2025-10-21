### simulation of rejection rate for different chainlengths

source("sim_rejection_rate.R")

#setup
A <- matrix(c(-0.4, 0, -0.3, 0.8), nrow = 2, byrow = T)
B <- 10
L <- 10

#different chainlengths
T500 <- 500
n500 <- T500 / L

message("T: ", T500)
df_T500 <- sim_rej_rate(T500, n500, L, B, A, seq(0.01, 1, 0.01))
saveRDS(df_T500, file = "T500.rds")

### ------------------ ####
T2K <- 2000
n2K <- T2K / L

message("T: ", T2K)
df_T2K <- sim_rej_rate(T2K, n2K, L, B, A, seq(0.01, 1, 0.01))
saveRDS(df_T2K, file = "T2K.rds")

### ------------------ ####
T5K <- 5000
n5K <- T5K / L

message("T: ", T5K)
df_T5K <- sim_rej_rate(T5K, n5K, L, B, A, seq(0.01, 1, 0.01))
saveRDS(df_T5K, file = "T5K.rds")

### ------------------ ####
#note: only ran with repetitions = 200 and not 500
T10K <- 10000
n10K <- T10K / L

message("T: ", T10K)
df_T10K <- sim_rej_rate(T10K, n10K, L, B, A, seq(0.01, 1, 0.01), repetitions = 200)
saveRDS(df_T10K, file = "T10K.rds")

gather_plots(df_T500)
gather_plots(df_T2K)
gather_plots(df_T5K)
gather_plots(df_T10K)


