#################### Hierarchical NS EVD ##########################################################################
# This code is for fitting a hierarchical non-stationary EVD model with all three variables maxWS, damages, minCP as 
# univariate EVD models.

### logpost_factory:
#This function computes the log unnormalized posterior density (log-likelihood + log-prior)
### Input: 
# priorscale: Scaling for hyperparameters
# data: (x2, x1, z1) i.e. (damages, maxWS, minCP) and covariates (avgLat, avgLon, Start month, Year, maximum category of storm)

### Metrop function in mcmc package:
## Use to run the MH-algorithm
## Input:
# step_scale: step scales for MH-steps
# pos_init: the initial state of the Markov chain

## Author: Sakshi Arya (aryax010@umn.edu)
## Date : 4th August, 2021, edited Jan 10, 2023 edited March 24, 2024
####################################################################################################################
rm(list = ls())
library(evd)
library(extRemes)
library(mcmc)
library(mcmcse)
library(xtable)
library(dplyr)



## File storage destinations
filep <- "/Users/sakshiarya/Library/CloudStorage/GoogleDrive-arya.sakshi44@gmail.com/My\ Drive/Research/2022/Ansu/BayesStorm/Tropical_Hurricane_Damages_Bayesian/2024_retake/"
filefigs <- "/Users/sakshiarya/Library/CloudStorage/GoogleDrive-arya.sakshi44@gmail.com/My\ Drive/Research/2022/Ansu/BayesStorm/Tropical_Hurricane_Damages_Bayesian/2024_retake/R\ Code\ Individual/Predictions/2020Prediction"
file.derived.21 <- "/Users/sakshiarya/Library/CloudStorage/GoogleDrive-arya.sakshi44@gmail.com/My\ Drive/Research/2022/Ansu/BayesStorm/Tropical_Hurricane_Damages_Bayesian/2024_retake/Data\ 2022"
file.Rcode <- "/Users/sakshiarya/Library/CloudStorage/GoogleDrive-arya.sakshi44@gmail.com/My\ Drive/Research/2022/Ansu/BayesStorm/Tropical_Hurricane_Damages_Bayesian/2024_retake/R\ Code\ Individual"

source(file.path(file.Rcode,"logpost_factory_allcovs.R"))
## Reading the data:
storms<-read.csv(file.path(file.derived.21,'Categorized_Storm_1960_2022.csv'))
data_cov <- read.delim(file.path(file.derived.21,'Annual_Covariates_1960_2022.csv'), header = TRUE, sep = ",")
head(data_cov)
data_cov <- data_cov[,-which(colnames(data_cov) == "YR")]



names(storms)
storm_damage<-storms[storms$minCP< Inf & storms$CURRENT.DAMAGE.2022>0 & !is.na(storms$CURRENT.DAMAGE.2022),]
head(storm_damage)
tail(storm_damage)

data_combined <- left_join(storm_damage, data_cov, by = "Year")
head(data_combined)
dim(data_combined)

data1960to2019 <- data_combined[data_combined$Year < 2020,]
dim(data1960to2019)
tail(data1960to2019)

tri_var<- data1960to2019[,c('minCP','maxWS','CURRENT.DAMAGE.2022')]
tri_var$logdamage<-log(tri_var$CURRENT.DAMAGE.2022)
tri_var$logmaxWS<-log(tri_var$maxWS)

tri_var$logminCP<-log(1013-tri_var$minCP)
tri_var$avgLat <- data1960to2019$avgLat
tri_var$avgLon <- data1960to2019$avgLon
tri_var$StartMonth <- data1960to2019$StartMonth
tri_var$Year <- data1960to2019$Year
tri_var$max_Cat <- data1960to2019$max_Cat
dim(tri_var)
tri_var$x1 <- tri_var$logmaxWS
tri_var$x2 <- tri_var$logdamage
tri_var$z1 <- tri_var$logminCP
tri_var[15:21] <- data1960to2019[,17:23]





# normalize covariates
head(tri_var)
tri_var <- tri_var[,-c(1:6)]
dim(tri_var)
tri_var <- tri_var[, c(1:5,9:15,6:8)]
tri_var[,1:12] <- scale(tri_var[,1:12])
head(tri_var)
apply(tri_var, 2, 'sd')
tri_var$scalex1 <- scale(tri_var$x1)
tri_var$scalez1 <- scale(tri_var$z1)

head(tri_var)
apply(tri_var,2, class)
tri_var_sub <- tri_var[,-which(colnames(tri_var) %in% c("max_Cat", "total"))]
head(tri_var_sub)

## Selecting the covariates we want in the model
tri_var <- tri_var_sub
covs1 <- colnames(tri_var[1:10])
d1 <- length(covs1)

## For the whole model
m_z1 <- fevd(z1, data = tri_var,location.fun=as.formula(paste("~",paste(covs1, collapse = "+"))), scale.fun=~1, shape.fun=~1,type='GEV')
m_z1$results$par
m_x1 <- fevd(x1, data = tri_var,location.fun=as.formula(paste("~scalez1 +", paste(covs1, collapse = "+"))), scale.fun=~1, shape.fun=~1,type='GEV')
m_x1$results$par
m_x2 <- fevd(x2, data = tri_var,location.fun=as.formula(paste("~ scalez1 + scalex1 +", paste(covs1, collapse = "+"))), scale.fun=~1, shape.fun=~1,type='GEV')
m_x2$results$par

freq_se <- c(summary(m_z1)$se.theta, summary(m_x1)$se.theta, summary(m_x2)$se.theta)
names(freq_se)

pos_init <- c(m_z1$results$par, m_x1$results$par, m_x2$results$par)
xtable(as.data.frame(pos_init), digits = 4)


pos_init <- as.numeric(c(pos_init))
is.numeric(pos_init)
# pos_init[42] <- -pos_init[42]

pos_init1 <- pos_init[c(1:(d1+1),(d1+4):((d1+4)+(d1+1)),(2*d1+5 + 3):(2*d1+5 + 3 +(d1+2)),(d1+3):(d1+2),(2*d1+5+2):(2*d1 + 6),(3*d1 + 12):(3*d1 + 11))]
pos_init1[c(37,39)] <- -0.4


names(pos_init1) <- c("Intercept minCP (Z1)", covs1[1:d1], 
                      "Intercept xaxWS (X1)", colnames(tri_var)[c(d1+5,1:d1)], 
                      "Intercept Damages (X2)",colnames(tri_var)[c(d1+5,d1+4,1:d1)],
                      expression(xi[Z1]),  expression(sigma[Z1]),
                      expression(xi[X1]),  expression(sigma[X1]),
                      expression(xi[X2]),  expression(sigma[X2]))
length(pos_init1)

## Scaling for hyperparameters
var_alpha <- rep(100, (d1+1));
alpha_sigz1 <- 1; beta_sigz1 <- 3; a_z1 <- -1; b_z1 <- 1
var_beta <- rep(1000, (d1+2)); alpha_sigx1 <- 1; beta_sigx1 <- 3;  a_x1 <- -0.55; b_x1 <- 0.5
var_gam <- c(rep(1000,2),rep(100,(d1+1))); alpha_sigx2 <- 1; beta_sigx2 <- 3; a_x2 <- -0.55; b_x2 <- 0.5


priorscale <- c(var_alpha,  
                var_beta, 
                var_gam, 
                alpha_sigz1, beta_sigz1, 
                alpha_sigx1, beta_sigx1, 
                alpha_sigx2, beta_sigx2,
                a_z1, b_z1,a_x1, b_x1, a_x2, b_x2)

length(priorscale)

step_scale <- c(rep(0.01, d1+1), rep(0.01, d1+2), rep(0.01, d1+3), 0.005,0.006,0.005, 0.006, 0.005, 0.006)
length(step_scale)
step_scale[c(25,26)] <- 0.01
# step_scale <- rep(0.05, length(pos_init1))
# step_scale[c(28:29,38)] <- c(0.012,0.01,0.05)
covs_model <- colnames(tri_var)[c(1:d1,(d1+4),(d1+5))]
length(covs_model)
logpost <- logpost_factory(data = tri_var, priorscale =  priorscale, covs = covs_model)
logpost(pos_init1)

set.seed(2403)
mo1 <- metrop(logpost, initial = pos_init1, 1e6, scale = step_scale)
mo1$accept
results <- mo1$batch

set.seed(2303)
TotalMC <- 1e6
t <- 1
reps <- 1000
results <- matrix(NA, ncol = length(pos_init1), nrow = TotalMC)
fac <- as.factor(seq(1,reps,by = 1))
indices <- rep(fac, each = TotalMC/reps)
groups <- split(1:TotalMC, indices)
while(t <= reps){
  mo1 <- metrop(logpost, initial = pos_init1, TotalMC/reps, scale = step_scale)
  results[c(groups[[t]]),] <- mo1$batch
  print(mo1$accept)
  if(mo1$accept < 0.20){
    step_scale <- 0.95*step_scale
  }else{
    step_scale <- 1.05*step_scale
  }
  pos_init1 <- mo1$batch[TotalMC/reps,]
  t <- t + 1
}

results <- mo1$batch
colnames(results) <- names(pos_init1)
save(results, file = file.path(filefigs,"hierNSallvariables_1e6_1960to2019.RData"))




load(file.path(filefigs,"hierNSallvariables_1e6_1960to2019.RData"))
post_means <- apply(results, 2, mean)
post_sd <- apply(results, 2, sd)

final_posterior_est <- cbind(post_means, post_sd)
names_row <- paste(names(pos_init1), c("", rep("MinCP",10), "", rep("MaxWS",11), "", rep("Damages",12), rep("MinCP",2), rep("MaxWS",2), rep("Damages",2)))
rownames(final_posterior_est) <- names_row
xtable(final_posterior_est, digits = 4)
xtable(cbind(post_means, post_sd, names(pos_init1)), digits = c(0,4,4,0))
results_every100 <- results[seq(1, nrow(results), 100),]
post_means_every100 <- apply(results_every100,2, mean)
post_sd_every100 <- apply(results_every100,2, sd)
xtable(cbind(post_means_every100, post_sd_every100, pos_init1), digits = 4)

dim(results_every100)
colnames(results_every100) <- names_row
cbind(post_means, pos_init1, post_sd)
cbind(post_means_every100, pos_init1, post_sd_every100)

pdf(file.path(filefigs,"TS_Hier_18April2023_all_1to6.pdf"))
plot.ts(results_every100[,1:6], main = "Time series plots for the MCMC output", cex.lab=0.6)
dev.off()
pdf(file.path(filefigs,"TS_Hier_18April2023_all_7to13.pdf"))
plot.ts(results_every100[,7:13], main = "Time series plots for the MCMC output", cex.lab=0.5)
dev.off()
pdf(file.path(filefigs,"TS_Hier_18April2023_all_14to21.pdf"))
plot.ts(results_every100[,14:21], main = "Time series plots for the MCMC output", cex.lab=0.6)
dev.off()
pdf(file.path(filefigs,"TS_Hier_18April2023_all_22to29.pdf"))
plot.ts(results_every100[,22:29], main = "Time series plots for the MCMC output", cex.lab=0.5)
dev.off()
pdf(file.path(filefigs,"TS_Hier_18April2023_all_30to36.pdf"))
plot.ts(results_every100[,30:36], main = "Time series plots for the MCMC output", cex.lab=0.5)
dev.off()
pdf(file.path(filefigs,"TS_Hier_18April2023_all_37to42.pdf"))
plot.ts(results_every100[,37:42], main = "Time series plots for the MCMC output", cex.lab = 0.6)
dev.off()

# 
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_1to5.pdf"))
acf(results_every100[,1:5], cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_6to10.pdf"))
acf(results_every100[,6:10], cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_11to15.pdf"))
acf(results_every100[,11:15], cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_16to20.pdf"))
acf(results_every100[,16:20], cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_21to25.pdf"))
acf(results_every100[,21:25], cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_26to30.pdf"))
acf(results_every100[,26:30], cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_31to35.pdf"))
acf(results_every100[,31:35], cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_36to40.pdf"))
acf(results_every100[,36:40],cex.main = 0.6)
dev.off()
pdf(file.path(filefigs,"ACF_Hier_18April2023_all_41to42.pdf"))
acf(results_every100[,41:42],cex.main = 0.6)
dev.off()

### Finding important variables:
depth_var <- rep(NA, length(pos_init1))
for(i in 1:ncol(results)){
  if(all(results[,i] >0) || all(results[,i] <0)){
    depth_var[i] <- 0
  }
  else{
    depth_var[i] <- 4*(table((results[,i]< 0))[2]/nrow(results))*(1-(table((results[,i]< 0))[2]/nrow(results)))
  }
}
depth_var

depth_var <- as.data.frame(depth_var)
names_row <- paste(names(pos_init1), c("", rep("MinCP",10), "", rep("MaxWS",11), "", rep("Damages",12), rep("MinCP",2), rep("MaxWS",2), rep("Damages",2)))
rownames(depth_var) <- names_row

depth_table <- xtable(data.frame(depth_var))
autoformat(depth_table)
display(depth_table)[3] <- "f"
digits(depth_table) <- c(0,4)
depth_table

digits(depth_table) <- xdigits(depth_table)
save(depth_table, file = file.path(filefigs,"DepthTablehierNSallvariables_1e6_1960to2019.RData"))

load(file.path(filefigs,"DepthTablehierNSallvariables_1e6_1960to2019.RData"))
depth_table




#############################################################################################
#########################################################################################
############## Repeat only using the selected covariates #############
## Selecting the covariates we want in the model, run until tri_var from before
tri_var <- tri_var_sub
tri_var_sub1 <- tri_var[1:10]
covs_minCP <- colnames(tri_var_sub1[,-which(colnames(tri_var_sub1) %in% c("NAO","SOI", "ANOM.3.4"))])
covs_damages <- colnames(tri_var_sub1[,-which(colnames(tri_var_sub1) %in% c("StartMonth","avgLon","SOI", "Year","NAO","Sunspots"))])


## For the whole model
m_z1 <- fevd(z1, data = tri_var,location.fun=as.formula(paste("~",paste(covs_minCP, collapse = "+"))), scale.fun=~1, shape.fun=~1,type='GEV')
m_z1$results$par
m_x1 <- fevd(x1, data = tri_var,location.fun=as.formula("~ scalez1"), scale.fun=~1, shape.fun=~1,type='GEV')
m_x1$results$par
m_x2 <- fevd(x2, data = tri_var,location.fun=as.formula(paste("~ scalez1 + scalex1 +", paste(covs_damages, collapse = "+"))), scale.fun=~1, shape.fun=~1,type='GEV')
m_x2$results$par

freq_se <- c(summary(m_z1)$se.theta, summary(m_x1)$se.theta, summary(m_x2)$se.theta)
names(freq_se)

pos_init <- c(m_z1$results$par, m_x1$results$par, m_x2$results$par)
xtable(as.data.frame(pos_init), digits = 4)


pos_init <- as.numeric(c(pos_init))
is.numeric(pos_init)
# pos_init[42] <- -pos_init[42]

pos_init1 <- pos_init[c(1:(length(covs_minCP)+1), (length(covs_minCP)+4):(length(covs_minCP)+5), (length(covs_minCP)+8): (length(covs_minCP)+ length(covs_damages)+10),
                        (length(covs_minCP)+3):(length(covs_minCP)+2), (length(covs_minCP)+7):(length(covs_minCP)+6), (length(covs_minCP)+length(covs_damages)+12):(length(covs_minCP)+length(covs_damages)+11))]


names(pos_init1) <- c("Intercept minCP (Z1)", covs_minCP, 
                      "Intercept maxWS (X1)", "scalez1",
                      "Intercept Damages (X2)", "scalex1", "scalez1",covs_damages,
                      expression(xi[Z1]),  expression(sigma[Z1]),
                      expression(xi[X1]),  expression(sigma[X1]),
                      expression(xi[X2]),  expression(sigma[X2]))
length(pos_init1)

## Scaling for hyperparameters
var_alpha <- rep(100, (length(covs_minCP)+1));
alpha_sigz1 <- 1; beta_sigz1 <- 3; a_z1 <- -1; b_z1 <- 1
var_beta <- rep(1000, 2); alpha_sigx1 <- 1; beta_sigx1 <- 3;  a_x1 <- -0.55; b_x1 <- 0.5
var_gam <- c(rep(1000,1),rep(100,(length(covs_damages)+2))); alpha_sigx2 <- 1; beta_sigx2 <- 3; a_x2 <- -0.55; b_x2 <- 0.5


priorscale <- c(var_alpha,  
                var_beta, 
                var_gam, 
                alpha_sigz1, beta_sigz1, 
                alpha_sigx1, beta_sigx1, 
                alpha_sigx2, beta_sigx2,
                a_z1, b_z1,a_x1, b_x1, a_x2, b_x2)

length(priorscale)

step_scale <- c(rep(0.01, length(covs_minCP)+1), rep(0.01,2), rep(0.01, length(covs_damages)+3), 0.005,0.006,0.005, 0.006, 0.005, 0.006)
length(step_scale)
# step_scale <- rep(0.05, length(pos_init1))
# step_scale[c(28:29,38)] <- c(0.012,0.01,0.05)


covs_model <- list(covs_minCP,"scalez1", c("scalez1","scalex1", covs_damages))

source(file.path(file.Rcode,"logpost_factory_selected_covs.R"))
logpost <- logpost_factory(data = tri_var, priorscale =  priorscale, covs = covs_model)
logpost(pos_init1)

set.seed(1706)
mo1 <- metrop(logpost, initial = pos_init1, 1e6, scale = step_scale)
mo1$accept
results <- mo1$batch


save(results, file = file.path(filefigs,"hierNSselecteddat_1e6_1960to2019.RData"))
colnames(results) <- names(pos_init1)
results <- mo1$batch

### To make the acceptance rate ~20%, we do the following
set.seed(2303)
TotalMC <- 1e6
t <- 1
reps <- 1000
results <- matrix(NA, ncol = length(pos_init1), nrow = TotalMC)
fac <- as.factor(seq(1,reps,by = 1))
indices <- rep(fac, each = TotalMC/reps)
groups <- split(1:TotalMC, indices)
while(t <= reps){
  mo1 <- metrop(logpost, initial = pos_init1, TotalMC/reps, scale = step_scale)
  results[c(groups[[t]]),] <- mo1$batch
  print(mo1$accept)
  if(mo1$accept < 0.23){
    step_scale <- 0.95*step_scale
  }else{
    step_scale <- 1.05*step_scale
  }
  pos_init1 <- mo1$batch[TotalMC/reps,]
  t <- t + 1
}

colnames(results) <- names(pos_init1)
head(results)
save(results, file = file.path(filefigs,"hierNSfulldat_1e6_23March2024.RData"))



