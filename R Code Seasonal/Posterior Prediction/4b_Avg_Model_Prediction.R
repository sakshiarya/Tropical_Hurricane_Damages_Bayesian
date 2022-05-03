#################################################################
#################################################################
#Title: Posterior Predictives using Avg. Estimate for 2016 based on 
#               1960-2015 using Method 2
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Create models and run MCMC
#Last Updated: 2/26/16
#Last updated: 4/6/21 (by Sakshi Arya)
#################################################################
#################################################################
#Clear out any R junk
#rm(list=ls())

#Loads all necessary libraries
library(mcmc)
library(mcmcse)
library(mvtnorm)
library(MCMCpack)
library(ggplot2)
library(grid)
library(reshape2)
library(parallel)
library(plyr)
library(dplyr)
library(tseries)
library(stringr)
library(dclone)
library(R2jags)
library(coda)
library(pscl)

#Run all data creation and hyperparameter code
filep<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Posterior Prediction/Predicting 2016/'
#source(paste(filep,'0_Data_Processing.r',sep=''))
#source(paste(filep,'1a_Hyperpar_Beta_2010.r',sep=''))
#source(paste(filep,'1b_Hyperpar_Theta_2010.r',sep=''))
#source(paste(filep,'1c_Hyperpar_Mu_Sigma_2010.r',sep=''))
#source(paste(filep,'2_Models_2010.r',sep=''))

filedata<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Data Sets/'
file_simulated_data <- '/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Posterior Prediction/Predicting 2016/Datasets/'
#Annualized data
Annual_Data<- read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Annual_1960_2019.csv'))
testsettruth<-Annual_Data[Annual_Data$Year>=2016,]

#Annualized Covariates (May/June Averages)
Annual_Cov1 <- read.csv(file.path(filedata,'Derived Data Sets 2021/Annual_Covariates_1960_2019.csv'))
Annual_Cov <- Annual_Cov1[Annual_Cov1$Year>=2016,]

#Observed covariate data for 2016
testset<- as.numeric(Annual_Cov[Annual_Cov$Year == 2016,4:9])

#Loading in the MCMC Chains
load(paste(filep,'MCMC Chains/clones.MCMC.final.FB.all.2015.rda',sep=''))
load(paste(filep,'MCMC Chains/jags.MCMC.final.FB.all.2015.rda',sep=''))

jags.2015.posterior <- as.mcmc(x= jags.MCMC.final.FB.all.2015)

#Saving the hierarchical estimates (means of the posterior)
pars2dclone<-summary(clones.MCMC.final.FB.all.2015)[[1]][,1]
pars2jags<-summary(jags.2015.posterior)$statistics[c(1:12,14:20),1]

#Assigning parameters for Group 1 (TS-2)
parameters_group1dclone<-c(pars2dclone['Beta1[1]'],pars2dclone['Beta1[2]'],pars2dclone['Beta1[3]'],pars2dclone['Beta1[4]'],pars2dclone['Beta1[5]'],pars2dclone['Beta1[6]'], pars2dclone['r'], pars2dclone['Theta1'], pars2dclone['log_mu1'], pars2dclone['log_sigma21'])
parameters_group1jags<-c(pars2jags['Beta1[1]'],pars2jags['Beta1[2]'],pars2jags['Beta1[3]'],pars2jags['Beta1[4]'],pars2jags['Beta1[5]'],pars2jags['Beta1[6]'],pars2jags['r'], pars2jags['Theta1'], pars2jags['log_mu1'], pars2jags['log_sigma21'])

#Assigning parameters for Group 2 (3-5)
parameters_group2dclone <-c(pars2dclone['Beta2[1]'], pars2dclone['Beta2[2]'],pars2dclone['Beta2[3]'],pars2dclone['Beta2[4]'],pars2dclone['Beta2[5]'],pars2dclone['Beta2[6]'], pars2dclone['Theta2'], pars2dclone['log_mu2'], pars2dclone['log_sigma22'])
parameters_group2jags <-c(pars2jags['Beta2[1]'], pars2jags['Beta2[2]'],pars2jags['Beta2[3]'],pars2jags['Beta2[4]'],pars2jags['Beta2[5]'],pars2jags['Beta2[6]'], pars2jags['Theta2'], pars2jags['log_mu2'], pars2jags['log_sigma22'])


posterior.sample.size<-dim(jags.MCMC.final.FB.all.2015$BUGSoutput$sims.list$Beta1)[1]
posterior.sample<-jags.MCMC.final.FB.all.2015$BUGSoutput$sims.list

# No. of simulations per MCMC sample 
num_simulations<-1000

#################################################################
#################################################################
#Category TS-2
#################################################################
#################################################################

#Making matrices to store simulated posterior predictive
N1<-matrix(nrow= num_simulations,ncol= posterior.sample.size)
L1<-matrix(nrow= num_simulations,ncol= posterior.sample.size)
D1<-matrix(nrow= num_simulations,ncol= posterior.sample.size)
lambdas<-matrix(nrow= num_simulations,ncol= posterior.sample.size)

#Simulating the posterior predictive for the Group 1 storms
set.seed(13238354)
start_time <- Sys.time()
for(j in 1:posterior.sample.size){
 lambdas[,j] <- exp(testset %*% posterior.sample$Beta1[j,])
 N1[,j]<- rnegbin(n=num_simulations, mu=lambdas[,j], theta=posterior.sample$r[j,1])
 L1[,j]<-rbinom(num_simulations,N1[,j], posterior.sample$Theta1[j,1])
 for(i in 1:num_simulations){
  if(L1[i,j]>0){
    D1[i,j]<-rnorm(1, posterior.sample$log_mu1[j,1],sqrt(posterior.sample$log_sigma21[j,1]))	
  } else{
    D1[i,j]<-0
  }
 }
}
end_time <- Sys.time()
end_time-start_time # about 10 mins

file_simulated_data <- '/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Posterior Prediction/Predicting 2016/Datasets/'
simulated_data_TS2_2015_4b <- list(N1=N1,L1=L1,D1=D1)
#save(simulated_data_TS2_2015_4b,file=file.path(file_simulated_data,'simulated_data_TS2_2015_4b.rda'))

#################################################################
#Frequency Calculations (TS-2)
#################################################################
#load(paste(file_simulated_data,"simulated_data_TS2_2015_4b.rda",sep=''))

## Frequency of occurence
tablemerge_each_N1 <- data.frame(matrix(NA,nrow = 100, ncol = 1))
freq_func <- function(i){
  tablemerge_each_N1[as.numeric(rownames(table(N1[,i])/sum(table(N1[,i]))))+1,2] <- table(N1[,i])/sum(table(N1[,i])) 
  return(tablemerge_each_N1)
}
start_time <- Sys.time()
result_N1 <- lapply(1:posterior.sample.size, freq_func)
end_time <- Sys.time()
end_time-start_time ## About 10 mins
result_Means_N1 <- rowMeans(do.call(cbind, result_N1), na.rm = TRUE)
#Posterior predictive density estimate for occurence:
result_Mean_N1 <- cbind(0:99, result_Means_N1)

#save(result_N1, file = "result_N1_mat_TS2.rda")
#save(result_Mean_N1, file = "N1_TS2_2016_meanposteriorpredict.rda")


## Landfall:
tablemerge_each_L1 <- data.frame(matrix(NA,nrow = 100, ncol = 1))
freq_func <- function(i){
  tablemerge_each_L1[as.numeric(rownames(table(L1[,i])/sum(table(L1[,i]))))+1,2] <- table(L1[,i])/sum(table(L1[,i])) 
  return(tablemerge_each_L1)
}
result_L1 <- lapply(1:posterior.sample.size, freq_func)
result_Means_L1 <- rowMeans(do.call(cbind, result_L1), na.rm = TRUE)
#Posterior predictive density estimate for landfall:
result_Mean_L1 <- cbind(0:99, result_Means_L1)


#save(result_L1, file = "result_L1_mat_TS2.rda")
#save(result_Mean_L1, file = "L1_TS2_2016_meanposteriorpredict.rda")


## Damage posterior predictive distribution:
avg.zeros<-mean(rowMeans(D1==0))
values<-data.frame(val=density(D1[,3][D1[,3]>0],bw=0.25,from=8, to=30,n=150)$x)
val.set_D1 <-data.frame(values,
                  prob=density(D1[, 1][D1[, 1]>0],bw=0.25,from=8, to=30,n=150)$y)

pdf_func <- function(i){
  return(density(D1[, i][D1[, i]>0],bw=0.25,from=8, to=30,n=150)$y)
}
result_D1 <- lapply(1:posterior.sample.size, pdf_func)
result_Means_D1 <- rowMeans(do.call(cbind, result_D1), na.rm = TRUE)
# Posterior predictive density for damages:
result_Mean_D1 <- cbind(values, result_Means_D1)


#2016 Simulated Frequency for TS-2 with observed value
load("Datasets/N1_TS2_2016_meanposteriorpredict.rda")
#pdf("Pred2016_Frequency_TS2_methodb.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(result_Mean_N1, type = "h",xlim=c(0,60),ylim=c(0,0.11),xlab='TS-2 Storm Frequency',main='',ylab='Density', cex.lab = 1.25, cex.axis = 1.25)
ts2.freq.truth.2016<-testsettruth[testsettruth$Year==2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq
segments(x0= ts2.freq.truth.2016 +0.1,y0=0, x1= ts2.freq.truth.2016 +0.1,y1=result_Mean_N1[which(result_Mean_N1[,1]== ts2.freq.truth.2016),2],col='red',lwd=3,lty=2)
segments(x0= ts2.freq.truth.2016-0.1,y0=0, x1= ts2.freq.truth.2016-0.1,y1=result_Mean_N1[which(result_Mean_N1[,1]== ts2.freq.truth.2016),2],col='red',lwd=3,lty=2)
#dev.off()


#################################################################
#Landfall Calculations (TS-2)
#################################################################
#2016 Simulated Landfall for TS-2 with observed value
load("Datasets/L1_TS2_2016_meanposteriorpredict.rda")
#pdf("Pred2016_Landfall_TS2_methodb.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(result_Mean_L1, type = "h",xlim=c(0,15), ylim = c(0,.34),cex.lab = 1.25,cex.axis = 1.25,xlab='TS-2 Landfall Frequency',main='',ylab='Density')
ts2.land.truth.2016<-testsettruth[testsettruth$Year==2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall
segments(x0= ts2.land.truth.2016 +0.1,y0=0, x1= ts2.land.truth.2016 +0.1,y1=result_Mean_L1[which(result_Mean_L1[,1]== ts2.land.truth.2016),2],col='red',lwd=3,lty=2)
segments(x0= ts2.land.truth.2016-0.1,y0=0, x1= ts2.land.truth.2016-0.1,y1=result_Mean_L1[which(result_Mean_L1[,1]== ts2.land.truth.2016),2],col='red',lwd=3,lty=2)
#dev.off()
#################################################################
#Damage Calculations (TS-2)
#################################################################
#2016 Simulated Damages for TS-2 with observed value
load("Datasets/D1_TS2_2016_meanposteriorpredict.rda")
load("Datasets/simulated_data_TS2_2015_4b.rda")
avg.zeros<-mean(rowMeans(simulated_data_TS2_2015_4b$D1==0))
#pdf("Pred2016_Damage_TS2_methodb.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(cbind(c(result_Mean_D1[,1],0),c((1-avg.zeros)*result_Mean_D1[,2], avg.zeros)),type='h',xlab='TS-2 Storm Damage',main='',ylab='Density', cex.lab = 1.25, cex.axis = 1.25, lwd = 1.5)
text(5,0.15,paste(round(avg.zeros,2),'chance of\n $0 Damage'), cex = 1.15)
text(25.4,0.17,paste(round(1-avg.zeros,2),'chance of\n log(Damage)\n in this distribution'), cex = 1.12)
ts2.damage.truth.2016 <- log(testsettruth[testsettruth$Year==2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage)
diff.ts2.damage.from.truth.2016 <- abs(ts2.damage.truth.2016 - result_Mean_D1[,1])
segments(x0= ts2.damage.truth.2016-0.1,y0=0, x1= ts2.damage.truth.2016-0.1,y1= result_Mean_D1[which(diff.ts2.damage.from.truth.2016 == min(diff.ts2.damage.from.truth.2016)),2],col='red',lwd=3,lty=2)
segments(x0= ts2.damage.truth.2016+0.1,y0=0, x1= ts2.damage.truth.2016+0.1,y1= result_Mean_D1[which(diff.ts2.damage.from.truth.2016 == min(diff.ts2.damage.from.truth.2016)),2],col='red',lwd=3,lty=2)
#dev.off()

#################################################################
#################################################################
#Category 3-5
#################################################################
#################################################################
#Making matrices to store simulated posterior predictive
N2<-matrix(nrow= num_simulations,ncol= posterior.sample.size)
L2<-matrix(nrow= num_simulations,ncol= posterior.sample.size)
D2<-matrix(nrow= num_simulations,ncol= posterior.sample.size)

#Simulating the posterior predictive for the Group 2 storms
set.seed(194370132)
for(j in 1:posterior.sample.size){
  N2[,j]<-rpois(num_simulations,exp(testset %*% posterior.sample$Beta2[j,]))
  L2[,j]<-rbinom(num_simulations,N2[,j], posterior.sample$Theta2[j])
 for(i in 1:num_simulations){
  if(L2[i,j]>0){
    D2[i,j]<-rnorm(1, posterior.sample$log_mu2[j],sqrt(posterior.sample$log_sigma22[j]))	
  } else{
    D2[i,j]<-0
  }
 }
}

simulated_data_35_2015_4b <- list(N2=N2,L2=L2,D2=D2)
#save(simulated_data_35_2015_4b,file=file.path(file_simulated_data,'simulated_data_35_2015_4b.rda'))

#################################################################
#Frequency Calculations (3-5)
#################################################################

## Frequency of occurence
tablemerge_each_N2 <- data.frame(matrix(NA,nrow = 100, ncol = 1))
freq_func <- function(i){
  tablemerge_each_N2[as.numeric(rownames(table(N2[,i])/sum(table(N2[,i]))))+1,2] <- table(N2[,i])/sum(table(N2[,i])) 
  return(tablemerge_each_N2)
}
result_N2 <- lapply(1:posterior.sample.size, freq_func)
result_Means_N2 <- rowMeans(do.call(cbind, result_N2), na.rm = TRUE)
#Posterior predictive density:
result_Mean_N2 <- cbind(0:99, result_Means_N2)
#save(result_N2, file = "result_N2_mat_35.rda")
#save(result_Mean_N2, file = "N2_35_2016_meanposteriorpredict.rda")

## Landfall:
tablemerge_each_L2 <- data.frame(matrix(NA,nrow = 100, ncol = 1))
freq_func_L2 <- function(i){
  tablemerge_each_L2[as.numeric(rownames(table(L2[,i])/sum(table(L2[,i]))))+1,2] <- table(L2[,i])/sum(table(L2[,i])) 
  return(tablemerge_each_L2)
}
#start_time <- Sys.time()
result_L2 <- lapply(1:posterior.sample.size, freq_func_L2)
result_Means_L2 <- rowMeans(do.call(cbind, result_L2), na.rm = TRUE)
#Posterior predictive density:
result_Mean_L2 <- cbind(0:99, result_Means_L2)
#end_time <- Sys.time()
#end_time-start_time # approx 5 min
#save(result_L2, file = "result_L2_mat_35.rda")
#save(result_Mean_L2, file = "L2_35_2016_meanposteriorpredict.rda")


## Damage posterior predictive distribution:
avg.zeros<-mean(rowMeans(D2==0))
values<-data.frame(val=density(D2[,3][D2[,3]>0],bw=0.25,from=8, to=30,n=150)$x)
val.set_D2 <-data.frame(values,
                        prob=density(D2[, 1][D2[, 1]>0],bw=0.25,from=8, to=30,n=150)$y)

pdf_func <- function(i){
  return(density(D2[, i][D2[, i]>0],bw=0.25,from=8, to=30,n=150)$y)
}
result_D2 <- lapply(1:posterior.sample.size, pdf_func)
result_Means_D2 <- rowMeans(do.call(cbind, result_D2), na.rm = TRUE)
#Posterior predictive density:
result_Mean_D2 <- cbind(values, result_Means_D2)

#save(result_D2, file = "result_D2_mat_35.rda")
#save(result_Mean_D2, file = "D2_35_2016_meanposteriorpredict.rda")

#2016 Simulated Frequency for 3-5 with observed value
#load("Datasets/N2_35_2016_meanposteriorpredict.rda")
#pdf("Pred2016_Frequency_T35_methodb.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(result_Mean_N2, type = "h",xlim=c(0,15), ylim = c(0,0.35),cex.lab = 1.25,cex.axis = 1.25, lwd = 1.5,xlab='3-5 Storm Frequency',main='',ylab='Density')
t35.freq.avgs <-testsettruth[testsettruth$Year==2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq
segments(x0= t35.freq.avgs +0.1,y0=0, x1= t35.freq.avgs +0.1,y1=result_Mean_N2[which(result_Mean_N2[,1]== t35.freq.avgs),2],col='red',lwd=3,lty=2)
segments(x0= t35.freq.avgs-0.1,y0=0, x1= t35.freq.avgs-0.1,y1=result_Mean_N2[which(result_Mean_N2[,1]== t35.freq.avgs),2],col='red',lwd=3,lty=2)
#dev.off()
#################################################################
#Landfall Calculations (3-5)
#################################################################

#2016 Simulated Landfall for 3-5 with observed value
load("Datasets/L2_35_2016_meanposteriorpredict.rda")
#pdf("Pred2016_Landfall_T35_methodb.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(result_Mean_L2, type = "h",xlim=c(0,10),cex.axis = 1.25,cex.lab = 1.25,lwd = 1.5,xlab='3-5 Landfall Frequency',main='',ylab='Density')
t35.land.truth.2016<-testsettruth[testsettruth$Year==2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall
segments(x0= t35.land.truth.2016 +0.1,y0=0, x1= t35.land.truth.2016 +0.1,y1=result_Mean_L2[which(result_Mean_L2[,1]== t35.land.truth.2016),2],col='red',lwd=3,lty=2)
segments(x0= t35.land.truth.2016-0.1,y0=0, x1= t35.land.truth.2016-0.1,y1=result_Mean_L2[which(result_Mean_L2[,1]== t35.land.truth.2016),2],col='red',lwd=3,lty=2)
#dev.off()
#################################################################
#Damage Calculations (3-5)
#################################################################

#2016 Simulated Damages for 3-5 with observed value
load("Datasets/D2_35_2016_meanposteriorpredict.rda")
load("Datasets/simulated_data_35_2015_4b.rda")
avg.zeros<-mean(rowMeans(simulated_data_35_2015_4b$D2==0))
#pdf("Pred2016_Damage_T35_methodb.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(cbind(c(result_Mean_D2[,1],0),c((1-avg.zeros)*result_Mean_D2[,2], avg.zeros)),type='h',xlab='3-5 Storm Damage',cex.lab = 1.25,cex.axis = 1.25,lwd = 1.5,main='',ylab='Average Posterior Predictive Density')
text(4.7,0.2,paste(round(avg.zeros,2),'chance of\n $0 Damage'), cex = 1.12)
text(25,0.20,paste(round(1-avg.zeros,2),'chance of\n log(Damage)\n in this distribution'), cex= 1.12)
t35.damage.truth.2016<-testsettruth[testsettruth$Year==2016 & testsettruth$Cat_HURDAT=='3-5',]$damage
diff.35.damage.from.truth.2016 <- abs(t35.damage.truth.2016 - result_Mean_D2[,1])
segments(x0= log(t35.damage.truth.2016)+0.1,y0=0, x1= log(t35.damage.truth.2016)+0.1,y1=0.057 ,col='red',lwd=3,lty=5)
segments(x0= log(t35.damage.truth.2016)-0.1,y0=0, x1= log(t35.damage.truth.2016)-0.1,y1=0.057 ,col='red',lwd=3,lty=5)
#dev.off()

##### Calculating Mahanalobis distance of the new observation from the PPD:####
# 2016 TS2:
load(paste(file_simulated_data,"simulated_data_TS2_2015_4b.rda",sep = ""))
dim(simulated_data_TS2_2015_4b$N1)
mean_vec <- sapply(simulated_data_TS2_2015_4b, mean)
colmean_vec <- sapply(simulated_data_TS2_2015_4b, rowMeans)
cov_NLD <- cov((colmean_vec))
mahal_dist <- mahalanobis(x=Annual_Data[Annual_Data$Year<2016,3:5],center= mean_vec, 
                          cov= cov_NLD)
depth <- 1/(1+mahal_dist)
mahal_dist_2016 <- mahalanobis(x=testsettruth[testsettruth$Year==2016&testsettruth$Cat_HURDAT=="TS-2",3:5],center= mean_vec, 
                               cov= cov_NLD)
# MB depth:
depth_2016 <- 1/(1+mahal_dist_2016)
#p-value (rank):
p_val_depth_TS2_2016 <- sum(depth < depth_2016)/length(depth)


########### T-35 2016 ###############
file_simulated_data <- '/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Posterior Prediction/Predicting 2016/Datasets/'
load(paste(file_simulated_data,"simulated_data_35_2015_4b.rda", sep = ""))
mean_vec2 <- sapply(simulated_data_35_2015_4b, mean)
colmean_vec2 <- sapply(simulated_data_35_2015_4b, rowMeans)
cov_NLD2 <- cov((colmean_vec2))
mahal_dist_T35 <- mahalanobis(x=Annual_Data[Annual_Data$Year<2016,3:5],
                              center= mean_vec2, cov= cov_NLD2)
depth_T35 <- 1/(1+mahal_dist_T35)
mahal_dist_2016_T35 <- mahalanobis(x=testsettruth[testsettruth$Year==2016&testsettruth$Cat_HURDAT=="3-5",3:5],center= mean_vec2, 
                               cov= cov_NLD2)
# MB depth for 2016 data point:
depth_2016_T35 <- 1/(1+mahal_dist_2016_T35)

#p-value/rank:
p_val_depth_T35_2016 <- sum(depth_T35 < depth_2016_T35)/length(depth_T35)
