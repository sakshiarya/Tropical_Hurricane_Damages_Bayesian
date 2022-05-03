#################################################################
#################################################################
#Title: Posterior Predictives using Best Estimate for 2016 based on 1960-2015
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Create models and run MCMC
#Created: 2/29/16
#Last updated: 4/6/21 (Sakshi Arya)
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
library(DepthProc)

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

#Observed Data from 2011 forward
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
num_simulations<-10^5

#Simulating the posterior predictive for the Group 1 storms
set.seed(99399)
 for(i in 1:num_simulations){
  lambdas[i] <- exp(testset %*% parameters_group1jags[1:6])
  #p_NB[i] <- parameters_group1jags['r']/(parameters_group1jags['r'] + lambdas)
  N1[i]<- rnegbin(n=1, mu=lambdas[i], theta=parameters_group1jags['r'])
  L1[i]<-rbinom(1,N1[i], parameters_group1jags['Theta1'])
  if(L1[i]>0){
    D1[i]<-rnorm(1, parameters_group1jags['log_mu1'],sqrt(parameters_group1jags['log_sigma21']))	
  } else{
    D1[i]<-0
  }
 }


#2016 Simulated Frequency for TS-2 with observed value
pdf("Pred2016_Frequency_TS2_methoda.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(table(N1)/sum(table(N1)), ylim = c(0,0.12), cex.lab = 1.25, cex.axis = 1.25, type = "h",ylab='Density',main='',xlab='TS-2 Storm Frequency',xlim=c(0,60),axes=F, lwd=1.5)
axis(side=2)
axis(side=1,at=seq(0,60,10))
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq +0.2,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq +0.2,y1= 0.085,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq-0.2,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Freq -0.2,y1= 0.085,col='red',lwd=3,lty=2)
dev.off()

#2016 Simulated Landfall for TS-2 with observed value
pdf("Pred2016_Landfall_TS2_methoda.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))
plot(table(L1)/sum(table(L1)), xlim=c(0,15), ylim = c(0,.34), cex.axis = 1.25, cex.lab = 1.25, lwd = 1.5, type = "h",ylab='Density',xlab='TS-2 Landfall Frequency',main='',axes=F)
axis(side=2, cex.axis = 1.25)
axis(side=1,at=seq(0,15,5),pos=-0.014,cex.axis=1.25)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall +0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall +0.1,y1= 0.335,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall -0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$Landfall -0.1,y1= 0.335,col='red',lwd=3,lty=2)
dev.off()

#2011 Simulated Damages for TS-2 with observed value
pdf("Pred2016_Damages_TS2_methoda.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(D1,freq=F,xlab='TS-2 Storm Damage',main='', breaks=c(0,1,2:31),axes=F,
     ylim=c(-0.01,0.32),cex.lab=1.25)
axis(side=2,at=seq(0,0.30,0.05),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,30,5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,31,0.5),pos=0.32,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.32,0.04),pos=31,lwd.ticks=0,labels=F)
text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage +0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage+0.1),y1= 0.115,col='red',lwd=3,lty=2)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage -0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage-0.1),y1= 0.115,col='red',lwd=3,lty=2)
segments(x0=0,y0=0.30,x1=0,y1=0.32)
dev.off()

#################################################################
#################################################################
#Category 3-5
#################################################################
#################################################################

#Making vectors to store simulated posterior predictive
N2 <-vector()
L2 <-vector()
D2 <-vector()
num_simulations<-10^5

#Simulating the posterior predictive for the Group 2 storms
set.seed(9919123)
for(i in 1:num_simulations){
  N2[i]<-rpois(1,exp(testset %*% parameters_group2jags[1:6]))
  L2[i]<-rbinom(1,N2[i], parameters_group2jags['Theta2'])
  if(L2[i]>0){
    D2[i]<-rnorm(1, parameters_group2jags['log_mu2'],sqrt(parameters_group2jags['log_sigma22']))	
  } else{
    D2[i]<-0
  }
}

#2016 Simulated Frequency for 3-5 with observed value
pdf("Pred2016_Frequency_T35_methoda.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(table(N2)/sum(table(N2)),cex.lab=1.25, cex.axis = 1.25, lwd = 1.5, type = "h",ylab='Density',main='',xlab='3-5 Storm Frequency',xlim=c(0,15), ylim = c(0,0.35),axes=F)
axis(side=2, cex.axis = 1.25)
axis(side=1,at=seq(0,15,5), cex.axis = 1.25)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq +0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq +0.1,y1= 0.042,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq-0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Freq -0.1,y1= 0.042,col='red',lwd=3,lty=2)
dev.off()

#2016 Simulated Landfall for 3-5 with observed value
pdf("Pred2016_Landfall_T35_methoda.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))
plot(table(L2)/sum(table(L2)),cex.lab=1.2, type = "h",ylab='Density',
     xlab='3-5 Landfall Frequency',main='',axes=F,xlim=c(0,10), cex.axis = 1.25)
axis(side=2, cex.axis = 1.25)
axis(side=1,at=seq(0,10,2), cex.axis = 1.25)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall +0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall +0.1,y1= 0.35,col='red',lwd=3,lty=2)
segments(x0= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall -0.1,y0=0, x1= testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$Landfall -0.1,y1= 0.35,col='red',lwd=3,lty=2)
dev.off()

#2016 Simulated Damages for 3-5 with observed value
pdf("Pred2016_Damages_T35_methoda.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(D2,freq=F,xlab='3-5 Storm Damage',main='',ylim=c(0,0.6), breaks=c(0,1,2:35),axes=F,cex.lab=1.25)
axis(side=2,pos=0)
axis(side=1,pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,35,1),pos=0.6,lwd.ticks=0,labels=F, cex.axis = 1.25)
axis(side=4,at=seq(0,0.6,0.01),pos=35,lwd.ticks=0,labels=F, cex.axis = 1.25)
text(8,0.3,paste(round(mean(D2==0),2),'chance of\n $0 Damage'),cex=1.2)
text(28,0.28,paste(round(1-mean(D2==0),2),'chance of\n log(Damage)\n in this distribution'),cex=1.2)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage+0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage+0.1),y1= 0.07,col='red',lwd=3,lty=2)
segments(x0= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage-0.1),y0=0, x1= log(testsettruth[testsettruth$Year== 2016 & testsettruth$Cat_HURDAT=='3-5',]$damage-0.1),y1= 0.07,col='red',lwd=3,lty=2)
dev.off()



# Data depth:
########################### Mahalanobis depth ####################################
##################################################################################
mahalanobisdepthrank <- function(dat, testdata){
  meanvec <- apply(dat, 2, mean)
  covNLD <- cov(dat)
  mahaldist <- mahalanobis(x=dat,center= meanvec, 
                           cov= covNLD)
  depth <- 1/(1+mahaldist)
  
  mahaldist_newdata <- mahalanobis(x= testdata,center= meanvec, 
                                   cov= covNLD)
  depth_newdata <- 1/(1+mahaldist_newdata)
  pvaldepth_newdata <- sum(depth < depth_newdata)/length(depth)
  return(c(depth_newdata, mahaldist_newdata, pvaldepth_newdata))
}



############################### Tukey Depth ######################################
##################################################################################

# Create own function:
Tukey_depth_rank <- function(dat = dat, reps, test_data = new_2016, u1 = u1){
mult <- matrix(NA, nrow = reps, ncol = dim(dat)[1])
prob_each <- matrix(NA, nrow = reps, ncol = dim(dat)[1])
for(i in 1:reps){
  mult[i,] <- c(u1[,i]%*%t(dat))
  prob_each[i,] <- (rank(mult[i,])-1)/dim(dat)[1]
}
depth_TS2 <- apply(prob_each, 2, min)
mult_test_data <- cbind(mult,  c(t(u1) %*% t(test_data)))
prob_each_mult_test_data <- matrix(NA, nrow = reps, ncol = dim(dat)[1]+ 1)
for(i in 1:reps){
  prob_each_mult_test_data[i,] <- (rank(mult_test_data[i,])-1)/(dim(dat)[1] +1)
}
depth_TS2_newtest_data <- min(prob_each_mult_test_data[,dim(dat)[1]+1])
## Compare using existing package:
depth_TS2_newtest_data_package <- depthTukey(as.numeric(test_data), dat)
pval_depth_test_data_TS2 <- sum(depth_TS2 < depth_TS2_newtest_data)/length(depth_TS2)
return(c(depth_TS2_newtest_data,depth_TS2_newtest_data_package,pval_depth_test_data_TS2))
}


# Generate observations on unit sphere
set.seed(189273)
reps <- 500
u <- matrix(rnorm(reps*3), ncol = 3)
norm_vec <- function(x) sqrt(sum(x^2))
u1 <- sapply(1:reps, function(i) {u[i,]/apply(u, 1, norm_vec)[i]})

# New data to find the depth for:
new_2016_TS2 <- testsettruth[testsettruth$Year==2016&testsettruth$Cat_HURDAT=="TS-2",3:5]
new_2016_TS2$damage <- log(new_2016_TS2$damage)
new_2016_T35 <- testsettruth[testsettruth$Year==2016&testsettruth$Cat_HURDAT=="3-5",3:5]
new_2016_T35$damage <- log(new_2016_T35$damage)


## Include 0 damages TS2 2016:
dat1 <- cbind(L1,N1, D1)
Tukey_depth_rank(dat = dat1, reps = reps, test_data = new_2016_TS2, u1 = u1)
mahalanobisdepthrank(dat = dat1, testdata = new_2016_TS2)

## Exclude 0 damges TS2:
dat2 <- cbind(L1[D1>0],N1[D1>0], D1[D1>0])
Tukey_depth_rank(dat = dat2, reps = reps, test_data = new_2016_TS2, u1 = u1)
mahalanobisdepthrank(dat = dat2, testdata = new_2016_TS2)

## Include 0 damages T35:
dat3 <- cbind(L2, N2, D2)
Tukey_depth_rank(dat = dat3, reps = reps, test_data = new_2016_T35, u1 = u1)
mahalanobisdepthrank(dat = dat3, testdata = new_2016_T35)

## Exclude 0 damages T35:
dat4 <- cbind(L2[D2>0],N2[D2>0], D2[D2>0])
Tukey_depth_rank(dat = dat4, reps = reps, test_data = new_2016_T35, u1 = u1)
mahalanobisdepthrank(dat = dat4, testdata = new_2016_T35)
