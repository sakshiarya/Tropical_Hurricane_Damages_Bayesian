#################################################################
#################################################################
#Title: Empirical Bayesian Hierarchical Models with Selected Covariates
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Create models and run MCMC
#Created: 2/25/16
#Last updated: 3/18/2021 (By Sakshi Arya)
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
library(rjags)
library(R2jags)
library(coda)
#install.packages("runjags")
library(runjags)
library(mvtnorm)
#library(jags)
#Run all data creation and hyperparameter code
filep<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/'
#source(paste(filep,'0_Data_Processing.r',sep=''))
source(file.path(filep,'1a_Hyperpar_Beta21.R'))
source(file.path(filep,'1b_Hyperpar_Theta21.R'))
source(file.path(filep,'1c_Hyperpar_Mu_Sigma21.r'))


#################################################################
#################################################################
#Data Sets used
#################################################################

filedata<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Data Sets/'
filemcmc<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/MCMC Chains/March 31/Empirical Bayesian/Selected_covs'

#Storm-wise data
Storm_Data <-read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Storm_1960_2019.csv'))

#Annualized data
Annual_Data<- read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Annual_1960_2019.csv'))

#Annualized Covariates (May/June Averages)
Annual_Cov <- read.csv(file.path(filedata,'Derived Data Sets 2021/Annual_Covariates_1960_2019.csv'))

Annual_Data <- Annual_Data[order(Annual_Data$Year),]

#################################################################
#################################################################
#BUGS model using May/June Avg. for the selected covariates
#Based on freq. model selection
#Uses NAO,AMO,ANOM 3.4 and Atl SST in TS-2, AMO, ANOM 3.4 and Atl_SST for 3-5
#Uses negative binomial specification for TS-2 hurricanes
#################################################################
#################################################################

model.bugs.onecov.nb <- function() {
  for (i in 1:num_obs) {
    num_storms1[i] ~ dnegbin(p[i],r)
    p[i] <- r/(r + exp(lambda1[i]))
    lambda1[i] <- inprod(Beta1[], model.mat.x1[i,])
    
    num_storms2[i] ~ dpois(exp(lambda2[i]))
    lambda2[i] <- inprod(Beta2[], model.mat.x2[i,])
    
    landfalling1[i] ~ dbin(Theta1, num_storms1[i])
    landfalling2[i] ~ dbin(Theta2, num_storms2[i])
    
    log_mu1_zero[i] <- (1-zero[i])*log_mu1 + zero[i]*1e-10
    log_tau1_zero[i] <- ((1-zero[i])^(-2))*log_tau1 
    damage_year1[i] ~ dnorm(log_mu1_zero[i],log_tau1_zero[i]) 
    damage_year2[i] ~ dnorm(log_mu2,log_tau2) 
    
    # zero inflation for Damage
    zero[i] ~ dbern((1-Theta1)^num_storms1[i])
}
  
  for(j in 1:4){
    Beta1[j] ~ dnorm(hyper_Beta1[1,j], hyper_Beta1[2,j])
  }
  for(j in 1:3){
    Beta2[j] ~ dnorm(hyper_Beta2[1,j], hyper_Beta2[2,j])
  }

  
  Theta1~ dbeta(hyper_Theta[1,1], hyper_Theta[1,2] )
  Theta2~ dbeta(hyper_Theta[2,1], hyper_Theta[2,2] )
  
  log_mu1 ~ dnorm(hyper_lognorm_Mu[1,1], hyper_lognorm_Mu[1,2])
  log_tau1 ~ dgamma(hyper_lognorm_Tau[1,1], hyper_lognorm_Tau[2,1])
  log_sigma21 <- 1/log_tau1
  
  log_mu2 ~ dnorm(hyper_lognorm_Mu[2,1], hyper_lognorm_Mu[2,2])
  log_tau2 ~ dgamma(hyper_lognorm_Tau[1,2], hyper_lognorm_Tau[2,2])
  log_sigma22 <- 1/log_tau2
  
  r ~ dunif(0,70)
}

#################################################################
#################################################################

#Putting in NA for 0 damages
Annual_Data[Annual_Data $Cat_HURDAT=='TS-2' & Annual_Data $Landfall==0,]$damage<-NA
Annual_Data[Annual_Data $Cat_HURDAT=='3-5' & Annual_Data $Landfall==0,]$damage<-NA
de <- data.frame(Year = 1962, Cat_HURDAT = "3-5", Landfall = 0, Freq = 0, damage = NA)
Annual_Data2 <- rbind(Annual_Data[Annual_Data$Cat_HURDAT=='3-5',], de)
Annual_Data2 <- Annual_Data2[order(Annual_Data2$Year),]

data.nbpoi.selected <-list(
  #Damages in data
  damage_year1 =  log(Annual_Data[Annual_Data$Cat_HURDAT=='TS-2' ,]$damage),
  damage_year2 =  log(Annual_Data2$damage), 
  #Covariates and storm frequency in data
  model.mat.x1= cbind(Annual_Cov[,c(4,6:8)]), 
  model.mat.x2= cbind(Annual_Cov[,c(6:8)]), 
  num_storms1= Annual_Data[Annual_Data$Cat_HURDAT=='TS-2' ,]$Freq, 
  num_storms2= Annual_Data2[Annual_Data2$Cat_HURDAT=='3-5',]$Freq, 
  
  #Landfalling frequency in data
  landfalling1= Annual_Data[Annual_Data$Cat_HURDAT=='TS-2' ,]$Landfall, 
  landfalling2= Annual_Data2[Annual_Data2$Cat_HURDAT=='3-5',]$Landfall, 
  
  #Prior specification for Beta/Phi (Frequency)
  hyper_Beta1 = hyper_Beta1_nb_selected[,-1], 
  hyper_Beta2 = hyper_Beta2_poi_selected[,-1], 
  
  #Prior specification for Theta (Landfall)
  hyper_Theta = hyper_Theta, #May need to change the [1,2] element to 1, otherwise node convergence issue
  
  #Prior specification for Mu/Tau (Damage)
  hyper_lognorm_Mu = hyper_lognorm_Mu,
  hyper_lognorm_Tau = hyper_lognorm_Tau, 
  
  #Total number of observations
  num_obs = dim(Annual_Cov)[1] )

#Parameters to be estimated
params.nb<-c("Beta1","r", "Beta2","Theta1", "Theta2", 'log_mu1','log_mu2', 'log_sigma21','log_sigma22')


#################################################################
#################################################################
#Model fitting dclone
#################################################################
#################################################################
n.clones=5;
n.chains=3;
n.iter= 1e5;
n.adapt =100;
n.thin=1;

set.seed(22415)
clones.MCMC.final.EB.select.April5.SST <-dc.fit(data= data.nbpoi.selected, params= params.nb, model= model.bugs.onecov.nb, n.iter = n.iter, n.clones= n.clones, multiply="num_obs", n.chains= n.chains,n.adapt = n.adapt, n.update = 0)
save(clones.MCMC.final.EB.select.April5.SST ,file=file.path(filemcmc,'clones.MCMC.final.EB.select.April5.SST.rda'))

#################################################################
#################################################################
#Model fitting JAGS
#################################################################
#################################################################
n.chains=3;
n.iter= 1e5;
n.burnin =100;
n.thin=1;

set.seed(22415)
jags.MCMC.final.EB.select.April5.SST <-jags(data= data.nbpoi.selected, parameters.to.save= params.nb, model.file= model.bugs.onecov.nb, n.chains= n.chains, n.thin= n.thin, n.iter= n.iter, n.burnin= n.burnin)
save(jags.MCMC.final.EB.select.April5.SST,file=file.path(filemcmc,'jags.MCMC.final.EB.select.April5.SST.rda'))



