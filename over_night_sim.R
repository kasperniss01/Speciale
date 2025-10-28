### simulation of rejection rate for different chainlengths
rm(list = ls())
source("sim_rejection_rate.R")
source("conditional_distributions.R")


#setup
set.seed(420)
A <- matrix(c(-0.4, -0.3, 0, 0.8), nrow = 2)
# A <- matrix(runif(4, -1, 1), nrow = 2)
# A[1, -1] <- 0
B <- 10
L <- 10


#different chainlengths
T500 <- 500

message("T: ", T500)
df_T500 <- sim_rej_rate(T500, L = 10, B = 10, A, seq(0.01, 1, 0.01),
                        parametric = TRUE,
                        remainder_true_ccfs = list(
                          true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                          true_psi = function(x, u, A, t) {
                            char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                          }
                        )
                        )
saveRDS(df_T500, file = "T500.rds")

### ------------------ ####
T1000 <- 1000

message("T: ", T1K)
df_T1K <- sim_rej_rate(T1000, L, B, A, seq(0.01, 1, 0.01),
                       parametric = TRUE, 
                       remainder_true_ccfs = list(
                         true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                         true_psi = function(x, u, A, t) {
                           char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                         }
                       ))
saveRDS(df_T1K, file = "T1K.rds")

### ------------------ ####
T2K <- 2000

message("T: ", T2K)
df_T2K <- sim_rej_rate(T2K, L, B, A, seq(0.01, 1, 0.01),
                       parametric = TRUE, 
                       remainder_true_ccfs = list(
                         true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                         true_psi = function(x, u, A, t) {
                           char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                         }
                       ))
saveRDS(df_T2K, file = "T2K.rds")

### ------------------ ####
T5K <- 5000
n5K <- T5K / L

message("T: ", T5K)
df_T5K <- sim_rej_rate(T5K, L, B, A, seq(0.01, 1, 0.01),
                       parametric = TRUE,
                       remainder_true_ccfs = list(
                         true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                         true_psi = function(x, u, A, t) {
                           char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                         }
                       ))
saveRDS(df_T5K, file = "T5K.rds")

### ------------------ ####
message("T: ", 7500)
df_T7500 <- sim_rej_rate(7500, L, B, A, seq(0.01, 1, 0.01),
                        parametric = TRUE,
                        remainder_true_ccfs = list(
                          true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                          true_psi = function(x, u, A, t) {
                            char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                          }
                        ))
saveRDS(df_T7500, file = "T7500.rds")

### ------------------ ####
T10K <- 10000

message("T: ", T10K)
df_T10K <- sim_rej_rate(T10K, L, B, A, seq(0.01, 1, 0.01),
                        parametric = TRUE,
                        remainder_true_ccfs = list(
                          true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                          true_psi = function(x, u, A, t) {
                            char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                          }
                        ))
saveRDS(df_T10K, file = "T10K.rds")

#rerunning usual stuff to also get remainders

### ----- higher dimensions ---- ###
#run code in high-dimension here
# not run this time

#Setup
set.seed(420)
A4D <- matrix(runif(16, min = -1, max = 1), nrow = 4, ncol = 4)
A4D[1, -1] <- 0

# set.seed(420)
# A2D <- matrix(runif(4, min = -1, max = 1), nrow = 2, ncol = 2)
# A2D[1, -1] <- 0
# 

#different chainlengths
T500 <- 500
now_running <- T500


message("T: ", T500, ", Y 4D")
df_T500_4D <- sim_rej_rate(500, L = 5, B = 10, A_matrix = A4D, alphas = seq(0.01, 1, 0.01),
                           remainder_true_ccfs = list(
                             true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                             true_psi = function(x, u, A, t) {char_func_cond_Y_given_X_highdim_mat(A, t, x, u)}
                             ),
                           parametric = TRUE
                        )

saveRDS(df_T500_4D, file = "T500_4D.rds")



### ------------------ ####
T2K <- 2000
now_running <- T2K

message("T: ", T2K , "Y 4D")
df_T2K_4D <- sim_rej_rate(T2K, L = 10, B = 10, A4D, seq(0.01, 1, 0.01),
                          remainder_true_ccfs = list(
                            true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                            true_psi = function(x, u, A, t) char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                            ),
                          parametric = TRUE
                          )

saveRDS(df_T2K_4D, file = "T2K_4D.rds")

### ------------------ ####
T5K <- 5000
now_running <- T5K

message("T: ", T5K, ", Y 4D")
df_T5K_4D <- sim_rej_rate(T5K, L = 10, B = 10, A4D, seq(0.01, 1, 0.01),
                          remainder_true_ccfs = list(
                            true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                            true_psi = function(x, u, A, t) char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                          ),
                          parametric = TRUE
                        )
saveRDS(df_T5K_4D, file = "T5K_4D.rds")

# ### ------------------ ####
T10K <- 10000
now_running <- T10K

message("T: ", T10K)
df_T10K_4D <- sim_rej_rate(T10K, L = 10, B = 10, A4D, seq(0.01, 1, 0.01),
                           remainder_true_ccfs = list(
                             true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                             true_psi = function(x, u, A, t) char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                           ),
                           parametric = TRUE
                          )
saveRDS(df_T10K_4D, file = "T10K_4D.rds")





