# ------------
# simulation 2: sensitivity analysis
# Methods: DR-PM, DR-ML
# ------------
rm(list = ls())
set.seed(123)
library(tidyverse)
library(geepack) #GEE
library(CRTgeeDR) #Aug-GEE
library(SuperLearner) # TMLE
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
mp <- 0.1 # Missing data proportion among covariates
pi <- 0.5 # Randomization probability
sim_size <- 1000

tictoc::tic()
results <- foreach(iter = 1:sim_size, .combine = cbind, .packages = package_list) %dopar% {
  # Complete data generation --------
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
    M <- sample(floor(N[j]/2)+(-2:2),1)
    cluster_id <- j
    data.frame(cbind(Y, R_Y, A=A[j], R_X1, X1, R_X2, X2, R_C1 = R_C1[j], C1=C1[j], N=N[j], M), cluster_id)
  })
  
  observed_data <- map_dfr(complete_data, function(d){
    obs_indi <- sample(1:d$N[1], size = d$M[1], replace = F)
    d[obs_indi,] %>% 
      mutate(RY = ifelse(R_Y, Y, NA), RX1 = ifelse(R_X1, X1, 0), RX2 = ifelse(R_X2, X2, 0), RC1 = ifelse(R_C1, C1, 0)) %>%
      dplyr::select(cluster_id, R_Y, RY, A, R_X1, RX1, R_X2, RX2, R_C1, RC1, M)
  })
  
  cl_data <- observed_data %>% group_by(cluster_id) %>% summarise_all(mean)
  
  tryCatch({
    # DR-PM -------
    #model treatment assignment
    pi_1 <- glm(A~ R_X1+RX1+R_X2+RX2+R_C1+RC1+M, family = binomial, data = cl_data) %>% predict(type = "response")
    pi_0 <- 1-pi_1
    # model missing outcomes
    gee_R_Y <- geeglm(R_Y~ A+R_X1+RX1+R_X2+RX2+R_C1+RC1+M, id = cluster_id, 
                      data =observed_data, family = "binomial", corstr = "independence")
    gee_R_Y_A1 <- predict(gee_R_Y, newdata = mutate(observed_data, A=1), type = "response")
    gee_R_Y_A0 <- predict(gee_R_Y, newdata = mutate(observed_data, A=0), type = "response")
    # outcome modeling
    gee_Y <- geeglm(RY~ A+R_X1+RX1+R_X2+RX2+R_C1+RC1+M, id = cluster_id, 
                    data =observed_data[observed_data$R_Y==1,], family = "gaussian", corstr = "exchangeable")
    gee_Y_A1 <- predict(gee_Y, newdata = mutate(observed_data, A=1), type = "response")
    gee_Y_A0 <- predict(gee_Y, newdata = mutate(observed_data, A=0), type = "response")
    # Compute the EIF
    U <- data.frame(cluster_id = observed_data$cluster_id,
                    U1 = observed_data$A/rep(pi_1, cl_data$M) * ifelse(observed_data$R_Y, observed_data$RY - gee_Y_A1, 0)/gee_R_Y_A1 + gee_Y_A1,
                    U0 = (1-observed_data$A)/rep(pi_0, cl_data$M) * ifelse(observed_data$R_Y, observed_data$RY - gee_Y_A0, 0)/gee_R_Y_A0 + gee_Y_A0) %>%
      group_by(cluster_id) %>%
      summarise_all(mean)
    DR_PM_est <- mean(U$U1 - U$U0)
    # Compute the variance correction because of nuisance parameter estimation for outcome missingness
    dd <- model.matrix(R_Y~A+R_X1+RX1+R_X2+RX2+R_C1+RC1+M, data = observed_data)
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
    DR_PM_1mRY <- 1- mean(gee_R_Y_A)
    
    # DR-ML --------
    # modeling missing outcomes
    ml_R_Y <- SuperLearner(observed_data$R_Y, X = observed_data[,c("A", "R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], 
                           family = "binomial", id=observed_data$cluster_id, SL.library = SLmethods, cvControl = list(V=5))
    ml_R_Y_A1 <- predict(ml_R_Y, newdata = mutate(observed_data[,c("A", "R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], A=1), type = "response")$pred
    ml_R_Y_A0 <- predict(ml_R_Y, newdata = mutate(observed_data[,c("A", "R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], A=0), type = "response")$pred
    # outcome modeling
    ml_Y_fitA1 <- SuperLearner(observed_data$RY[observed_data$R_Y==1&observed_data$A==1], 
                               X = observed_data[observed_data$R_Y==1&observed_data$A==1, c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], 
                               family = "gaussian", id=observed_data$cluster_id[observed_data$R_Y==1&observed_data$A==1], 
                               SL.library = SLmethods, cvControl = list(V=5))
    ml_Y_A1 <- predict(ml_Y_fitA1, newdata = observed_data[,c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], type = "response")$pred
    ml_Y_fitA0 <- SuperLearner(observed_data$RY[observed_data$R_Y==1&observed_data$A==0], 
                               X = observed_data[observed_data$R_Y==1&observed_data$A==0, c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], 
                               family = "gaussian", id=observed_data$cluster_id[observed_data$R_Y==1&observed_data$A==0], 
                               SL.library = SLmethods, cvControl = list(V=5))
    ml_Y_A0 <- predict(ml_Y_fitA0, newdata = observed_data[,c("R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], type = "response")$pred
    # EIF
    U <- data.frame(cluster_id = observed_data$cluster_id,
                    U1 = observed_data$A/pi * ifelse(observed_data$R_Y, observed_data$RY - ml_Y_A1, 0)/ml_R_Y_A1 + ml_Y_A1,
                    U0 = (1-observed_data$A)/(1-pi) * ifelse(observed_data$R_Y, observed_data$RY - ml_Y_A0, 0)/ml_R_Y_A0 + ml_Y_A0) %>%
      group_by(cluster_id) %>%
      summarise_all(mean)
    DR_ML_est <- mean(U$U1 - U$U0)
    DR_ML_se <- sd(U$U1 - U$U0)/sqrt(m-7)
    DR_ML_1mRY <- 1- mean(predict(ml_R_Y, newdata = observed_data[,c("A", "R_X1", "RX1", "R_X2", "RX2", "R_C1", "RC1", "M")], type = "response")$pred)
    
    c(DR_PM_est, DR_PM_se, DR_PM_1mRY, DR_ML_est, DR_ML_se, DR_ML_1mRY)
  }, error = function(e) { return(rep(NA, 6))})
}
tictoc::toc()
stopCluster(cl)

sensitivity <- expand.grid(m = m, delta1m0 = c(0,1,2,3,4), method = c("CDR-PM", "CDR-ML")) %>% mutate(mean = NA, sd = NA)
for(i in 1:nrow(sensitivity)){
  if(sensitivity$method[i] == "CDR-PM"){
    ci.lower <- results[1,] - qt(0.975, m) * results[2,]
    gamma1m0 <- (ci.lower - sensitivity$delta1m0[i]/2)/results[3,]
  } else {
    ci.lower <- results[4,] - qt(0.975, m) * results[5,]
    gamma1m0 <- (ci.lower - sensitivity$delta1m0[i]/2)/results[6,]
  }
  gamma1m0[abs(gamma1m0)>100] <- NA
  sensitivity$mean[i] <- mean(gamma1m0, na.rm =T)
  sensitivity$sd[i] <- sd(gamma1m0, na.rm=T)
}
sensitivity
saveRDS(sensitivity, file = paste0("sim2-m",m,"-mp",mp,".rds"))


# making plot------
plot_d <- readRDS("sim2-m30-mp0.1.rds")
p1 <- ggplot(plot_d) +
  geom_point(aes(y=mean, x=delta1m0, color=method), size=2, position=position_dodge(.3)) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, x = delta1m0, color=method), width=.2, position=position_dodge(.3)) +
  theme_bw() +
  labs(x = expression(delta[1]-delta[0]), y = expression(gamma[1]-gamma[0])) +
  ylim(-10,45)+
  ggtitle(paste0("m=30, 10% missing data"))+
  theme(text=element_text(size = 16,family = "Times New Roman"), legend.position = c(0.85, 0.85))


plot_d <- readRDS("sim2-m100-mp0.1.rds")
p2 <- ggplot(plot_d) +
  geom_point(aes(y=mean, x=delta1m0, color=method), size=2, position=position_dodge(.3)) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, x = delta1m0, color=method), width=.2, position=position_dodge(.3)) +
  theme_bw() +
  ylim(-10,45)+
  labs(x = expression(delta[1]-delta[0]), y = expression(gamma[1]-gamma[0])) +
  ggtitle(paste0("m=100, 10% missing data"))+
  theme(text=element_text(size = 16,family = "Times New Roman"), legend.position = c(0.85, 0.85))

plot_d <- readRDS("sim2-m30-mp0.3.rds")
p3<- ggplot(plot_d) +
  geom_point(aes(y=mean, x=delta1m0, color=method), size=2, position=position_dodge(.3)) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, x = delta1m0, color=method), width=.2, position=position_dodge(.3)) +
  theme_bw() +
  ylim(-5,17)+
  labs(x = expression(delta[1]-delta[0]), y = expression(gamma[1]-gamma[0])) +
  ggtitle(paste0("m=30, 30% missing data"))+
  theme(text=element_text(size = 16,family = "Times New Roman"), legend.position = c(0.85, 0.85))


plot_d <- readRDS("sim2-m100-mp0.3.rds")
p4 <- ggplot(plot_d) +
  geom_point(aes(y=mean, x=delta1m0, color=method), size=2, position=position_dodge(.3)) +
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd, x = delta1m0, color=method), width=.2, position=position_dodge(.3)) +
  theme_bw() +
  ylim(-5,17)+
  labs(x = expression(delta[1]-delta[0]), y = expression(gamma[1]-gamma[0])) +
  ggtitle(paste0("m=100, 30% missing data"))+
  theme(text=element_text(size = 16,family = "Times New Roman"), legend.position = c(0.85, 0.85))

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, labels = c('(a)', '(b)', '(c)','(d)'), label_size = 16)
cowplot::save_plot("sim-sensitivity.png",p, base_width = 10, base_height = 10)
