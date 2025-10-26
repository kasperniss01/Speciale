### simulation of rejection rate for different chainlengths

source("sim_rejection_rate.R")

#deprecated now?

#setup
A <- matrix(c(-0.4, 0, -0.3, 0.8), nrow = 2, byrow = T)
B <- 10
L <- 10

#different chainlengths
T500 <- 500
n500 <- T500 / L

message("T: ", T500)
df_T500 <- sim_rej_rate(T500, L, B, A, seq(0.01, 1, 0.01))
saveRDS(df_T500, file = "T500.rds")

### ------------------ ####
T2K <- 2000
n2K <- T2K / L

message("T: ", T2K)
df_T2K <- sim_rej_rate(T2K, L, B, A, seq(0.01, 1, 0.01))
saveRDS(df_T2K, file = "T2K.rds")

### ------------------ ####
T5K <- 5000
n5K <- T5K / L

message("T: ", T5K)
df_T5K <- sim_rej_rate(T5K, L, B, A, seq(0.01, 1, 0.01))
saveRDS(df_T5K, file = "T5K.rds")

### ------------------ ####
T10K <- 10000
n10K <- T10K / L

message("T: ", T10K)
df_T10K <- sim_rej_rate(T10K, L, B, A, seq(0.01, 1, 0.01))
saveRDS(df_T10K, file = "T10K.rds")

#rerunning usual stuff to also get remainders

### ----- higher dimensions ---- ###
#run code in high-dimension here
# not run this time

#Setup
A4D <- matrix(runif(16, min = -1, max = 1), nrow = 4, ncol = 4)
A4D[1, 2:4] <- 0

B <- B
L <- L

#different chainlengths
T500 <- 500

message("T: ", T500, ", Y 4D")
df_T500_4D <- sim_rej_rate(T500, L, B, A4D, seq(0.01, 1, 0.01))
saveRDS(df_T500_4D, file = "T500_4D.rds")

### ------------------ ####
T2K <- 2000

message("T: ", T2K , "Y 4D")
df_T2K_4D <- sim_rej_rate(T2K, L, B, A4D, seq(0.01, 1, 0.01))
saveRDS(df_T2K_4D, file = "T2K_4D.rds")

### ------------------ ####
T5K <- 5000

message("T: ", T5K, ", Y 4D")
df_T5K_4D <- sim_rej_rate(T5K, L, B, A4D, seq(0.01, 1, 0.01))
saveRDS(df_T5K_4D, file = "T5K_4D.rds")

# ### ------------------ ####
# T10K <- 10000
# 
# message("T: ", T10K)
# df_T10K_4D <- sim_rej_rate(T10K, L, B, A4D, seq(0.01, 1, 0.01))
# saveRDS(df_T10K, file = "T10K.rds")





