source("sim_rejection_rate.R")

### 4-dimensional (i.e. 3d Y)
alternative_matrix <- function(gamma, d = 3, vec = rep(1, d), seed  = NULL) {
  vec <- vec/norm(vec, type = "2")
  
  if(!is.null(seed)) set.seed(seed)
  
  A <- matrix(runif((d+1)^2, -1,1), nrow = d+1, ncol = d+1)
  
  A <- (A / abs(eigen(A)$values[1])) %>% round(digits = 2)
  
  
  A[1,-1] <- gamma*vec
  
  A  
}

gammas <- seq(0,1,0.333)
powerTs <- c(2000, 3000)

PowerAnalysisLists<- list()

for(Tlen in powerTs) {
  PowerAnalysisLists[[paste0("T",as.character(Tlen))]] <- list()
  
  
  for(gamma in gammas) {
    A_alt <- alternative_matrix(gamma = gamma, d = 3, vec = c(1,1,1), seed = 100)
    
    message("T: ", Tlen, " | gamma: ", gamma)
    df_power <- sim_rej_rate(Tlen, L = 10, B = 10, A = A_alt, alphas = seq(0.01, 1, 0.01),
                            parametric = TRUE,
                            remainder_true_ccfs = list(
                              true_phi = function(x, u, A) char_func_cond_X_next_given_X_previous_mat(A, x, u),
                              true_psi = function(x, u, A, t) {
                                char_func_cond_Y_given_X_highdim_mat(A, t, x, u)
                              }
                            )
                          )
    
    PowerAnalysisLists[[paste0("T",as.character(Tlen))]][[paste0("gamma", as.character(gamma))]] <- df_power
  }
  
}

saveRDS(PowerAnalysisLists, file = "datasets/KasperPowerAnalysis_AR1_4D.rds")




