# ------------
# simulation 1: continuous outcome, no sampling, MAR
# Methods: unadjusted, IPW, DR, DR-PM, DR-ML
# ------------
rm(list = ls())
set.seed(123)
library(tidyverse)
library(geepack) #GEE
library(CRTgeeDR) #Aug-GEE
library(SuperLearner) # model prediction
expit <- function(x){1/(1+exp(-x))}
logit <- function(x){log(x/(1-x))}

library(foreach)
library(doSNOW)
cl <- makeCluster(6)
registerDoSNOW(cl)
SLmethods <- c("SL.glm", "SL.gam", "SL.rpart")
package_list <- c("tidyverse","geepack", "CRTgeeDR", "SuperLearner")

# Parameters setting ------------
m <- 100 # Number of clusters
mp <- 0.3 # Missing data proportion among covariates
pi <- 0.5 # Randomization probability
sim_size <- 1000

tictoc::tic()
results <- foreach(iter = 1:sim_size, .combine = cbind, .packages = package_list) %dopar% {
  # Complete data generation
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
  
  tryCatch({
  
  # Unadjusted-----
  sy <- observed_data %>% group_by(cluster_id) %>% summarise(RY_bar = mean(RY, na.rm=T), A = mean(A))
  unadj_est <- mean(sy$RY_bar[sy$A==1]) - mean(sy$RY_bar[sy$A==0])
  unadj_se <- sqrt(var(sy$RY_bar[sy$A==1])/pi + var(sy$RY_bar[sy$A==0])/(1-pi))/sqrt(m-1)
  
  # IPW-----
  ipw <-geeDREstimation(formula=RY~A, nameY = "RY", nameMISS= "R_Y", nameTRT = "A",
                        id="cluster_id" , data = observed_data,
                        family = "gaussian", corstr = "independence",
                        sandwich.nuisance = T,
                        model.weights=I(R_Y==1)~A + R_X1 + RX1 + R_X2 + RX2 + R_C1 + RC1 + N)
  ipw_est <- summary(ipw)$beta[2]
  ipw_se <- summary(ipw)$se.robust[2]
  
  # DR-aug-GEE----
  aug_gee <- geeDREstimation(formula=RY~A, nameY = "RY", nameMISS= "R_Y", nameTRT = "A",
                             id="cluster_id" , data = observed_data,
                             family = "gaussian", corstr = "independence",
                             sandwich.nuisance = T, stepwise.augmentation=FALSE,
                             model.weights=I(R_Y==1)~A + R_X1 + RX1 + R_X2 + RX2 + R_C1 + RC1 + N,
                             model.augmentation.trt= RY~ R_X1 + RX1 + R_X2 + RX2 + R_C1 + RC1 + N,
                             model.augmentation.ctrl= RY~ R_X1 + RX1 + R_X2 + RX2 + R_C1 + RC1 + N,
                             fay.adjustment = T
  )
  aug_gee_est <- summary(aug_gee)$beta[2]
  aug_gee_se <- summary(aug_gee)$se.robust[2]
  
  # DR-PM -------
  #model treatment assignment
  pi_1 <- glm(A~ R_X1+RX1+R_X2+RX2+R_C1+RC1+N, family = binomial, data = cl_data) %>% predict(type = "response")
  pi_0 <- 1-pi_1
  # model missing outcomes
  gee_R_Y <- geeglm(R_Y~ A+R_X1+RX1+R_X2+RX2+R_C1+RC1+N, id = cluster_id, 
                    data =observed_data, family = "binomial", corstr = "independence")
  gee_R_Y_A1 <- predict(gee_R_Y, newdata = mutate(observed_data, A=1), type = "response")
  gee_R_Y_A0 <- predict(gee_R_Y, newdata = mutate(observed_data, A=0), type = "response")
  # outcome modeling
  gee_Y <- geeglm(RY~ A+R_X1+RX1+R_X2+RX2+R_C1+RC1+N, id = cluster_id, 
                  data =observed_data[observed_data$R_Y==1,], family = "gaussian", corstr = "exchangeable")
  gee_Y_A1 <- predict(gee_Y, newdata = mutate(observed_data, A=1), type = "response")
  gee_Y_A0 <- predict(gee_Y, newdata = mutate(observed_data, A=0), type = "response")
  # Compute the EIF
  U <- data.frame(cluster_id = observed_data$cluster_id,
                  U1 = observed_data$A/rep(pi_1, cl_data$N) * ifelse(observed_data$R_Y, observed_data$RY - gee_Y_A1, 0)/gee_R_Y_A1 + gee_Y_A1,
                  U0 = (1-observed_data$A)/rep(pi_0, cl_data$N) * ifelse(observed_data$R_Y, observed_data$RY - gee_Y_A0, 0)/gee_R_Y_A0 + gee_Y_A0) %>%
    group_by(cluster_id) %>%
    summarise_all(mean)
  DR_PM_est <- mean(U$U1 - U$U0)
  # Compute the variance correction because of nuisance parameter estimation for outcome missingness
  dd <- model.matrix(R_Y~A+R_X1+RX1+R_X2+RX2+R_C1+RC1+N, data = observed_data)
  gee_R_Y_A <- predict(gee_R_Y, type = "response")
  D <- map(cl_data$cluster_id, function(j){
    index <- which(observed_data$cluster_id ==j)
    t(dd[index,]) %*% diag(gee_R_Y_A[index] * (1-gee_R_Y_A[index])) %*% dd[index,]
  }) %>% reduce(., `+`)/m
  l1 <- map(cl_data$cluster_id, function(j){
    index <- which(observed_data$cluster_id ==j)
    tt <- diag( (1-gee_R_Y_A1[index])/gee_R_Y_A1[index] * ifelse(observed_data$R_Y[index], observed_data$RY[index] - gee_Y_A1[index], 0)) %*% dd[index,]
    colMeans(tt)
  }) %>% reduce(., `+`)/m
  l0 <- map(cl_data$cluster_id, function(j){
    index <- which(observed_data$cluster_id ==j)
    tt <- diag( (1-gee_R_Y_A0[index])/gee_R_Y_A0[index] * ifelse(observed_data$R_Y[index], observed_data$RY[index] - gee_Y_A0[index], 0)) %*% dd[index,]
    colMeans(tt)
  }) %>% reduce(., `+`)/m
  correction1 <- t(l1 - l0) %*% solve(D) %*% (l1 - l0)
  #  Compute the variance correction because of nuisance parameter estimation for treatment assignment
  dd <- model.matrix(A~R_X1+RX1+R_X2+RX2+R_C1+RC1+N, data = cl_data)
  D <- pi * (1-pi) * t(dd) %*% dd
  l1 <- t(aggregate(ifelse(observed_data$R_Y, observed_data$RY - gee_Y_A1, 0)/gee_R_Y_A1, by = list(id=observed_data$cluster_id), mean)[,2]) %*% dd *(1-pi)/m
  l0 <- t(aggregate(ifelse(observed_data$R_Y, observed_data$RY - gee_Y_A0, 0)/gee_R_Y_A0, by = list(id=observed_data$cluster_id), mean)[,2]) %*% dd *(-pi)/m
  correction2 <- (l1 - l0) %*% solve(D) %*% t(l1 - l0)
  # Output se estimate
  DR_PM_se <- sqrt(var(U$U1-U$U0) - correction1 - correction2)/sqrt(m-7)
  

  # DR-ML --------
  # modeling missing outcomes
  ml_R_Y <- SuperLearner(observed_data$R_Y, X = observed_data[,c("A", "R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "N")], 
                         family = "binomial", id=observed_data$cluster_id, SL.library = SLmethods, cvControl = list(V=5))
  ml_R_Y_A1 <- predict(ml_R_Y, newdata = mutate(observed_data[,c("A", "R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "N")], A=1), type = "response")$pred
  ml_R_Y_A0 <- predict(ml_R_Y, newdata = mutate(observed_data[,c("A", "R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "N")], A=0), type = "response")$pred
  # outcome modeling
  ml_Y_fitA1 <- SuperLearner(observed_data$RY[observed_data$R_Y==1&observed_data$A==1], 
                             X = observed_data[observed_data$R_Y==1&observed_data$A==1, c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "N")], 
                             family = "gaussian", id=observed_data$cluster_id[observed_data$R_Y==1&observed_data$A==1], 
                             SL.library = SLmethods, cvControl = list(V=5))
  ml_Y_A1 <- predict(ml_Y_fitA1, newdata = observed_data[,c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "N")], type = "response")$pred
  ml_Y_fitA0 <- SuperLearner(observed_data$RY[observed_data$R_Y==1&observed_data$A==0], 
                             X = observed_data[observed_data$R_Y==1&observed_data$A==0, c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "N")], 
                             family = "gaussian", id=observed_data$cluster_id[observed_data$R_Y==1&observed_data$A==0], 
                             SL.library = SLmethods, cvControl = list(V=5))
  ml_Y_A0 <- predict(ml_Y_fitA0, newdata = observed_data[,c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "N")], type = "response")$pred
  # EIF
  U <- data.frame(cluster_id = observed_data$cluster_id,
                  U1 = observed_data$A/pi * ifelse(observed_data$R_Y, observed_data$RY - ml_Y_A1, 0)/ml_R_Y_A1 + ml_Y_A1,
                  U0 = (1-observed_data$A)/(1-pi) * ifelse(observed_data$R_Y, observed_data$RY - ml_Y_A0, 0)/ml_R_Y_A0 + ml_Y_A0) %>%
    group_by(cluster_id) %>%
    summarise_all(mean)
  DR_ML_est <- mean(U$U1 - U$U0)
  DR_ML_se <- sd(U$U1 - U$U0)/sqrt(m-7)

  c(unadj_est, unadj_se, ipw_est, ipw_se, aug_gee_est, aug_gee_se, DR_PM_est, DR_PM_se, DR_ML_est, DR_ML_se)
  }, error = function(e) { return(rep(NA, 10))})
}
tictoc::toc()
stopCluster(cl)
saveRDS(results, paste0("simulation-1-1-m=",m,"-mp=",mp,".rds"))

# summarizing results
results <- readRDS(paste0("simulation-1-1-m=",m,"-mp=",mp,".rds"))
est_index <- c(1,3,5,7,9)
se_index <- c(2,4,6,8,10)
delta <- 5
results[7, abs(results[7,]) > 100] <- NA
results[8, abs(results[8,]) > 100] <- NA
results[6,] <- results[6,] * sqrt(m)/sqrt(m-7) # add finite-sample adjustment for DR-Aug-GEE

summary_results <- data.frame(methods = c("Unadj", "IPW", "DR", "CDR-PM", "CDR-ML"),
                              bias = apply(results[est_index,], 1, mean, na.rm = T) - delta,
                              emp_se = apply(results[est_index,], 1, sd, na.rm = T),
                              avg_se = apply(results[se_index,], 1, mean, na.rm = T),
                              cp = apply(abs((results[est_index,]-delta)/results[se_index,]) <= qt(0.975, m-7), 1, mean, na.rm = T)) %>% mutate_if(is.numeric, round,2)
summary_results
sum(is.na(results))
