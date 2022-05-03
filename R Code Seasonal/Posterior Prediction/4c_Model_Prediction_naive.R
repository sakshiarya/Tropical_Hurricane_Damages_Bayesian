#################################################################
#################################################################
#Title: Posterior Predictives using each MCMC iteration once for 2016 
#      based on 1960-2015 as in Method 3
#Author: Sakshi Arya (aryax010@umn.edu)
#Objective: Create models and run MCMC
# Last updated: 04/08/2021
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
library(plyr)
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
#################################################################
#################################################################
#Category TS-2
#################################################################
#################################################################

#Making vectors to store simulated posterior predictive
N1 <-vector()
L1 <-vector()
D1 <-vector()
lambdas <-vector()

#Simulating the posterior predictive for the Group 1 storms
set.seed(99399)
for(i in 1:posterior.sample.size){
  lambdas[i] <- exp(testset %*% posterior.sample$Beta1[i,])
  N1[i]<- rnegbin(1, mu=lambdas[i], theta=posterior.sample$r[i,1])
  L1[i]<-rbinom(1,N1[i], posterior.sample$Theta1[i,1])
  if(L1[i]>0){
    D1[i]<-rnorm(1, posterior.sample$log_mu1[i,1],sqrt(posterior.sample$log_sigma21[i,1]))	
  } else{
    D1[i]<-0
  }
}


#2016 Simulated Frequency for TS-2 with observed value
#pdf("Pred2016_Frequency_TS2_methodc.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(table(N1)/sum(table(N1)), type = "h",ylab='Density',main='',xlab='TS-2 Storm Frequency',xlim=c(0,60),axes=F, cex.lab = 1.25)
axis(side=2, pos = -2.5, cex.axis = 1.25)
axis(side=1,at=seq(0,60,10), cex.axis = 1.25, pos=-0.004)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq +0.2,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq +0.2,y1= 0.08,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq-0.2,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq -0.2,y1= 0.08,col='red',lwd=3,lty=2)
#dev.off()

#2016 Simulated Landfall for TS-2 with observed value
#pdf("Pred2016_Landfall_TS2_methodc.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(table(L1)/sum(table(L1)), ylim = c(0,0.35), type = "h",ylab='Density',xlab='TS-2 Landfall Frequency',main='',xlim=c(0,15),axes=F, cex.lab =1.25)
axis(side=2,pos = -0.6, cex.axis = 1.25)
axis(side=1,at=seq(0,15,5),cex.axis = 1.25, pos=-0.015)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall +0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall +0.1,y1= 0.33,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall -0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall -0.1,y1= 0.33,col='red',lwd=3,lty=2)
#dev.off()

#2016 Simulated Damages for TS-2 with observed value
#pdf("Pred2016_Damages_TS2_methodc.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(D1,freq=F,xlab='TS-2 Storm Damage',main='', breaks=c(0,1,2:31),axes=F,
     ylim=c(0,0.32),cex.lab=1.25)
axis(side=2,at=seq(0,0.32,0.05),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,31,5),pos=-0.001,cex.axis=1.25)
axis(side=3,at=seq(0,31,0.5),pos=0.31,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.31,0.01),pos=31,lwd.ticks=0,labels=F)
text(6,0.17,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
text(25,0.2,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage+0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage+0.1),y1= 0.11,col='red',lwd=3,lty=2)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage-0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage-0.1),y1= 0.11,col='red',lwd=3,lty=2)
segments(x0=0,y0=0.30,x1=0,y1=0.31)
#dev.off()
#################################################################
#################################################################
#Category 3-5
#################################################################
#################################################################

#Making vectors to store simulated posterior predictive
N2 <-vector()
L2 <-vector()
D2 <-vector()


#Simulating the posterior predictive for the Group 2 storms
set.seed(9919123)
for(i in 1:posterior.sample.size){
  N2[i]<-rpois(1,exp(testset %*% posterior.sample$Beta2[i,]))
  L2[i]<-rbinom(1,N2[i], posterior.sample$Theta2[i])
  if(L2[i]>0){
    D2[i]<-rnorm(1, posterior.sample$log_mu2[i],sqrt(posterior.sample$log_sigma22[i]))	
  } else{
    D2[i]<-0
  }
}

#2016 Simulated Frequency for 3-5 with observed value
#pdf("Pred2016_Frequency_T35_methodc.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(table(N2)/sum(table(N2)),cex.lab=1.25, type = "h",ylab='Density',main='',xlab='3-5 Storm Frequency',xlim=c(0,15),axes=F)
axis(side=2, cex.axis = 1.25)
axis(side=1,at=seq(0,15,5), cex.axis = 1.25)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq +0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq +0.1,y1= 0.048,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq-0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq -0.1,y1= 0.048,col='red',lwd=3,lty=2)
#dev.off()

#2016 Simulated Landfall for 3-5 with observed value
#pdf("Pred2016_Landfall_T35_methodc.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(table(L2)/sum(table(L2)),cex.lab=1.25, type = "h",ylab='Density',
     xlab='3-5 Landfall Frequency',main='',axes=F,xlim=c(0,10))
axis(side=2, cex.axis = 1.25)
axis(side=1,at=seq(0,10,2), cex.axis = 1.25)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall +0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall +0.1,y1= 0.32,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall -0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall -0.1,y1= 0.32,col='red',lwd=3,lty=2)
#dev.off()

#2016 Simulated Damages for 3-5 with observed value
#pdf("Pred2016_Damages_T35_methodc.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(D2,freq=F,xlab='3-5 Storm Damage',main='',ylim=c(0,0.6), axes = F, breaks=c(0,1,2:35),cex.lab=1.2)
axis(side=2,pos=0)
axis(side=1,pos=0,cex.axis=1.2)
axis(side=3,at=seq(0,35,1),pos=0.6,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.6,0.01),pos=35,lwd.ticks=0,labels=F)
text(8,0.3,paste(round(mean(D2==0),2),'chance of\n $0 Damage'),cex=1.2)
text(28,0.28,paste(round(1-mean(D2==0),2),'chance of\n log(Damage)\n in this distribution'),cex=1.2)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage+0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage+0.1),y1= 0.0455,col='red',lwd=3,lty=1)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage-0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage-0.1),y1= 0.0455,col='red',lwd=3,lty=1)
#segments(x0=0,y0=-0.01,x1=0,y1=0)
#segments(x0=35,y0=-0.01,x1=36,y1=-0.01)
#segments(x0=0,y0=0.30,x1=0,y1=0.32)
#dev.off()

par(mar=c(3,3,-0.1,-0.1)+0.1, mgp =c(2,1,0))#sets margins of plotting area
hist(D1,freq=F,xlab='TS-2 Storm Damage',main='', breaks=c(0,1,2:31),axes=F,
     ylim=c(-0.01,0.26),cex.lab=1.25)
axis(side=2,at=seq(0,0.25,0.05),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,31,5),pos=-0.01,cex.axis=1.25)
axis(side=3,at=seq(0,31,1),pos=0.26,lwd.ticks=0,labels=F)
axis(side=4,at=seq(-0.01,0.26,0.01),pos=31,lwd.ticks=0,labels=F)
text(6,0.12,paste(round(mean(D1==0),2),'chance of\n $0 Damage'))
text(26,0.115,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'))
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage+0.6,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage+0.6,y1= 0.25,col='red',lwd=3,lty=2)
segments(x0=0,y0=-0.01,x1=0,y1=0)
segments(x0=30,y0=-0.01,x1=31,y1=-0.01)
segments(x0=0,y0=0.25,x1=0,y1=0.26)

############################ Mahalanobis Data depth ###################################
## TS2
mean_vec_TS2 <- c(mean(N1),mean(L1),mean(D1))
cov_NLD_TS2 <- cov(cbind(N1,L1,D1))
mahal_dist_TS2 <- mahalanobis(x=Annual_Data[Annual_Data$Year<2016,3:5],center= mean_vec_TS2, 
                              cov= cov_NLD_TS2)
depth_TS2 <- 1/(1+mahal_dist_TS2)
mahal_dist_2016_TS2 <- mahalanobis(x=testsettruth[testsettruth$Year==2016&testsettruth$Cat_HURDAT=="TS-2",3:5],center= mean_vec_TS2, 
                                   cov= cov_NLD_TS2)
depth_2016_TS2 <- 1/(1+mahal_dist_2016_TS2)
pval_depth_2016_TS2 <- sum(depth_TS2 < depth_2016_TS2)/length(depth_TS2)

## T35
mean_vec_T35 <- c(mean(N2),mean(L2),mean(D2))
cov_NLD_T35 <- cov(cbind(N2,L2,D2))
mahal_dist_T35 <- mahalanobis(x=Annual_Data[Annual_Data$Year<2016,3:5],
                              center= mean_vec_T35, cov= cov_NLD_T35)
depth_T35 <- 1/(1+mahal_dist_T35)
mahal_dist_2016_T35 <- mahalanobis(x=testsettruth[testsettruth$Year==2016&testsettruth$Cat_HURDAT=="3-5",3:5],center= mean_vec_T35, 
                                   cov= cov_NLD_T35)
depth_2016_T35 <- 1/(1+mahal_dist_2016_T35)
pval_depth_2016_T35 <- sum(depth_T35 < depth_2016_T35)/length(depth_T35)

