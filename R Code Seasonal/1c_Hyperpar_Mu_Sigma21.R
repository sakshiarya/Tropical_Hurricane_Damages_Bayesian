#################################################################
#################################################################
#Title: Hyperparameters for Damages
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Determine hurricane damage hyperparameters, mu
#			and tau
#Created: 2/25/16
#Last edited: 3/19/2021 (by Sakshi)
#################################################################
#################################################################
#rm(list = ls())
#Loads all necessary libraries
library(lme4)
library(car)
library(ggplot2)

#################################################################
#################################################################
#Data Sets used
#################################################################
#################################################################

filedata<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Data Sets/'

#Storm-wise data
Storm_Data <-read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Storm_1960_2019.csv'))

#Annualized data
Annual_Data<- read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Annual_1960_2019.csv'))

#################################################################
#################################################################
#Hyperparameters for mu
#################################################################
#################################################################

#Annual Storm Damage >0
Annual_gt0_Damages<- Annual_Data[Annual_Data$damage>0,]
Annual_gt0_Damages$Cat_HURDAT <- as.factor(Annual_gt0_Damages$Cat_HURDAT)

Annual_gt0_Damages$Cat_HURDAT<-relevel(Annual_gt0_Damages$Cat_HURDAT,'TS-2')

#Initial Test of difference in hurricane categories 
summary(lm(log(damage)~Cat_HURDAT,data=Annual_gt0_Damages))

#Storm Wise Damage >0
Storm_Data1 <- Storm_Data[complete.cases(Storm_Data$CURRENT.DAMAGE.2021),]
Storm_gt0_Damages <- Storm_Data1[Storm_Data1$CURRENT.DAMAGE.2021 > 0,]
Storm_gt0_Damages$Cat_HURDAT <- as.factor(Storm_gt0_Damages$Cat_HURDAT)
Storm_gt0_Damages$Cat_HURDAT<-relevel(Storm_gt0_Damages$Cat_HURDAT,'TS-2')

#Using Annual Values to Create Hyperparameters
hyper_lognorm_Mu_mean<- with(Annual_gt0_Damages,tapply(damage, Cat_HURDAT,function(x){mean(log(x),na.rm=T)}))

#The stormwise values are not much different
tapply(Storm_Data1[Storm_Data1$CURRENT.DAMAGE.2021>0,]$CURRENT.DAMAGE.2021, Storm_Data1[Storm_Data1$CURRENT.DAMAGE.2021>0,]$Cat_HURDAT,function(x){mean(log(x),na.rm=T)})

#Setting the hyperparameters as the mean and then picking a 
#small value for 1/variance in the prior for mu
hyper_lognorm_Mu<-cbind(hyper_lognorm_Mu_mean, rep(1e-4,2))

#Density Plot of the Annual log Damages by Category
#pdf("Log_Damage_new.pdf")
ggplot(Annual_gt0_Damages,aes(x=log(damage),fill= Cat_HURDAT))+ geom_density(alpha=.25)+xlab('Annual log(Damage)')+ylab('Density')+ xlim(10, 30)+ scale_fill_manual('Max. Category',values=c("grey","darkblue"))+theme_bw()
#dev.off()

#Density Plot of the Stormwise log Damages by Category
ggplot(Storm_gt0_Damages,aes(x=log(CURRENT.DAMAGE.2021),fill= Cat_HURDAT))+ geom_density(alpha=.25)+xlab('Storm-wise log(Damage)')+ylab('Density')+ scale_fill_manual('Max. Category',values=c("grey","darkblue"))+theme_bw()

#Final priors for mu with hyperparameters specified
par(mar=c(3,3,0,0)+0.1, mgp =c(2,1,0))
curve(dnorm(x,hyper_lognorm_Mu[1,1],sqrt(1/hyper_lognorm_Mu[1,2])),from=-500,to=600,xlab= expression (paste(mu[1])),ylab='Density')
curve(dnorm(x,hyper_lognorm_Mu[2,1],sqrt(1/hyper_lognorm_Mu[2,2])),from=-500,to=600,xlab= expression (paste(mu[2])),ylab='Density')

#################################################################
#################################################################
#Hyperparameters for tau=1/sigma^2
#Setting the hyperparameters for Tau
#Shape is the sample Tau; Scale is set to 1
#This gives a mean and variance equal to the sample Tau
#################################################################
#################################################################

sample_tau_damages<-with(Annual_gt0_Damages,tapply(damage, Cat_HURDAT,function(x){1/var(log(x),na.rm=T)}))

#Setting the hyperparameters for Tau
#Shape is the sample Tau; Scale is set to 1
#This gives a mean and variance equal to the sample Tau
shape<- sample_tau_damages
scale<- rep(1,2)

hyper_lognorm_Tau<-rbind(shape, scale)

#Final priors for tau with hyperparameters specified
par(mar=c(3,3,0,0)+0.1, mgp =c(2,1,0))#sets margins of plotting area
curve(dgamma(x, shape= shape[1], scale= scale[1]),ylab='Density',from=0,to=5)
curve(dgamma(x, shape= shape[2], scale =scale[2]),ylab='Density',from=0,to=5)
