# ------------
# simulation 1: continuous outcome, no sampling, MAR
# Methods: two-stage tmle
# ------------
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)
param <- as.numeric(args[1])  # Takes the first argument as the iteration parameter
set.seed(param)
library(tidyverse)
library(SuperLearner) # model prediction
library(tmle) # tmle package
expit <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}
SLmethods <- c('SL.mean', 'SL.glm', 'SL.gam')

# Parameters setting ------------
p.grid <- expand.grid(m = c(30,100), mp = c(0.1, 0.3))
pi <- 0.5 # Randomization probability
results <- rep(NA, 8)

for(iter in 1: nrow(p.grid)){
  m <- p.grid[iter, 1] # Number of clusters
  mp <- p.grid[iter, 2] # Missing data proportion among covariates
  
  # data generation --------------
  N <- sample(10:90, size = m, replace = T)
  C1 <- rnorm(m, mean = N/10, sd = 1)
  R_C1 <- rbinom(m, size = 1, prob= expit(logit(1-mp) + 0.5*(C1 - N/10)))
  while(max(R_C1) == 0 | min(R_C1)>0.99){R_C1 <- rbinom(m, size = 1, prob= expit(logit(1-mp) + 0.5*(C1 - N/10)))}
  A <- rbinom(m, size = 1, prob = pi)
  complete_data <- map(1:m, function(j) {
    X1 <- rbinom(N[j], size = 1, prob = N[j]/100) 
    R_X1 <- rbinom(N[j], size = 1, prob= 1-mp)
    while(max(R_X1) == 0 | min(R_X1)>0.99){R_X1 <- rbinom(N[j], size = 1, prob= 1-mp)}
    X2 <- rnorm(N[j], mean = C1[j] * mean(X1), sd = 1) + rnorm(1, mean = 0, sd = 1) * (C1[j] > 0)
    R_X2 <- rbinom(N[j], size = 1, prob= expit(logit(1-mp) + 0.5*(C1[j] - N[j]/10)))
    while(max(R_X2) == 0| min(R_X2)>0.99){R_X2 <- rbinom(N[j], size = 1, prob= expit(logit(1-mp) + 0.5*(C1[j] - N[j]/10)))}
    Y <- sin(C1[j]) * exp(X1) * abs(1 + X2)/4 + A[j] * 10 * X1 + rnorm(1, mean = 0, sd = 1) +  rnorm(N[j], sd = 1)
    R_Y <- rbinom(N[j], size = 1, prob= expit(logit(0.99-0.2*mp)-(1.5+5*mp)*R_X1*X1))
    while(max(R_Y)==0| min(R_Y)>0.99){R_Y <- rbinom(N[j], size = 1, prob= expit(logit(0.99-0.2*mp)-(1.5+5*mp)*R_X1*X1))}
    # R_Y <- rbinom(N[j], size = 1, prob= expit(logit(0.95-1.5*mp) + 3 * R_X1*X1))
    # while(max(R_Y)==0| min(R_Y)>0.99){R_Y <- rbinom(N[j], size = 1, prob= expit(logit(0.95-1.5*mp) + 3 * R_X1*X1))}
    cluster_id <- j
    data.frame(cbind(Y, R_Y, A=A[j], R_X1, X1, R_X2, X2, R_C1 = R_C1[j], C1=C1[j], N=N[j]), cluster_id)
  })
  
  observed_data <- map_dfr(complete_data, function(d){
    d %>% 
      mutate(RY = ifelse(R_Y, Y, NA), RX1 = ifelse(R_X1, X1, 0), RX2 = ifelse(R_X2, X2, 0), RC1 = ifelse(R_C1, C1, 0)) %>%
      dplyr::select(cluster_id, R_Y, RY, A, R_X1, RX1, R_X2, RX2, R_C1, RC1, N)
  })
  
  cl_data <- observed_data %>% group_by(cluster_id) %>% summarise_all(mean)
  
  # run two-stage tmle --------------
  tryCatch({
    
    # Two-stage TMLE
    for(j in cl_data$cluster_id){
      Oj <- observed_data[observed_data$cluster_id == j, ]
      tmle.fit <- tmle(Y = Oj$RY, A = Oj$A, W = Oj[,c('R_X1', 'RX1', 'R_X2', 'RX2')], Delta = Oj$R_Y, 
                       Q.SL.library = SLmethods, V.Q=2,
                       g.SL.library = SLmethods, V.g=2,
                       g.Delta.SL.library = SLmethods, V.Delta=2)
      cl_data$RY[j] <- tmle.fit$estimates$EY1$psi
    }
    tmle.fit.cl <- tmle(Y = cl_data$RY, A = cl_data$A, W = cl_data[,c('R_X1', 'RX1', 'R_X2', 'RX2', "R_C1", "RC1", "N")], Q.SL.library = SLmethods, g.SL.library = SLmethods, V.Q=2, V.g=2)
    tmle_est <- tmle.fit.cl$estimates$ATE$psi
    tmle_se <- sqrt(tmle.fit.cl$estimates$ATE$var.psi)
    
    results[(2 * iter - 1):(2 * iter)] <- c(tmle_est, tmle_se)
  }, error = function(e) {results[(2 * iter - 1):(2 * iter)] <- return(rep(NA, 2))})
}



write.csv(results, file = paste0("temp-results/sim1-1-results_", param, ".csv"), row.names = FALSE)
# write.table(results, file = "results_combined.csv", append = TRUE, 
#             sep = ",", row.names = FALSE, col.names = (param == 1))

# run the following code to combine results -----------
files <- list.files(pattern = 'sim1-1-results_.*\\.csv')
results <- do.call(rbind, lapply(files, function(j) {t(read.csv(j))}))
write.csv(results, '../sim1-1-results.csv', row.names = FALSE)


# summarizing results ---------
results <- read.csv("TMLE/sim1-1-tmle-results.csv")
delta <- 5
est_index <- c(1,3,5,7)
se_index <- c(2,4,6,8)

summary_results <- data.frame(bias = apply(results[,est_index], 2, mean, na.rm = T) - delta,
                      emp_se = apply(results[,est_index], 2, sd, na.rm = T),
                      avg_se = apply(results[,se_index], 2, mean, na.rm = T), 
                      cp = NA)
summary_results$cp[c(1,3)] <- apply(abs((results[,c(1,5)]-delta)/results[,c(2,6)]) <= qt(0.975, 30-7), 2, mean, na.rm = T)
summary_results$cp[c(2,4)] <- apply(abs((results[,c(3,7)]-delta)/results[,c(4,8)]) <= qt(0.975, 100-7), 2, mean, na.rm = T)

round(summary_results, 2)

