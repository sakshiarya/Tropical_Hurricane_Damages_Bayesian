#################################################################
#################################################################
#Title: Posterior Predictives for Tri-variate EVD model using Best Estimate for 2016 based on 1960-2015
#Author: Sakshi Arya (aryax010@umn.edu)
#Objective: Create models and run MCMC
#Created: 08/27/21
#Last updated: 08/27/21
#################################################################
#################################################################
#Clear out any R junk
rm(list=ls())

#Loads all necessary libraries
library(mcmc)
library(mcmcse)
library(mvtnorm)
library(ggplot2)
library(grid)
library(reshape2)
library(plyr)
library(evd)
library(extRemes)

#Run all data creation and hyperparameter code
filep<-"/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/BayesianBivEVDMCMC/Predictions/Predicting 2016/Tri_EVD_predict2016/"
filedata<-"/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/BayesianBivEVDMCMC/Predictions/Predicting 2016"

#Annualized data
Annual_Data<- read.csv(file.path(filedata,'Categorized_Storm_1960_2019.csv'))
testsettruth<-Annual_Data[Annual_Data$Year==2016,]

storm_damage_2016 <- testsettruth[testsettruth$minCP< Inf & testsettruth$CURRENT.DAMAGE.2021>0 & !is.na(testsettruth$CURRENT.DAMAGE.2021),]

tri_var <-storm_damage_2016[,c('minCP','maxWS','CURRENT.DAMAGE.2021')]
tri_var$logdamage<-log(tri_var$CURRENT.DAMAGE.2021)
tri_var$logmaxWS<-log(tri_var$maxWS)

tri_var$logminCP<-log(1013-tri_var$minCP)
tri_var$avgLat <-storm_damage_2016$avgLat

#densityPlot(tri_var2015$logdamage)

tri_var$x1 <- tri_var$logmaxWS
tri_var$scalex1 <- scale(tri_var$logmaxWS)
tri_var$x2 <- log(tri_var$CURRENT.DAMAGE.2021)
tri_var$z1 <- tri_var$logminCP
tri_var$scalez1 <- scale(tri_var$logminCP)
tri_var$z2 <- scale(tri_var$avgLat)

# #Annualized Covariates (May/June Averages)
# Annual_Cov1 <- read.csv(file.path(filedata,'Derived Data Sets 2021/Annual_Covariates_1960_2019.csv'))
# Annual_Cov <- Annual_Cov1[Annual_Cov1$Year>=2016,]

#testset<- as.numeric(Annual_Cov[Annual_Cov$Year == 2016,4:9])

#Loading in the MCMC Chains
load(paste(filep,'NS_trivar_1e6_28Aug_1960to2015.RData',sep=''))



#Saving the hierarchical estimates (means of the posterior)
pars <- colMeans(results)


#Making vectors to store simulated posterior predictive
#lambdas <-vector()
num_simulations<-10^5
n <- dim(storm_damage_2016)[1]
vars_tri1 <-matrix(NA, ncol = 3, nrow = num_simulations)
vars_tri2 <-matrix(NA, ncol = 3, nrow = num_simulations)


#try
## list with marginals
marg_CP <- as.matrix(cbind(pars[1] + pars[2]*tri_var$z2, rep(pars[3],n),rep(pars[4], n)))
marg_WS <- as.matrix(cbind(pars[5] + pars[6]*tri_var$scalez1 + pars[7]*tri_var$z2, 
                           rep(pars[8],n), rep(pars[9],n)))
marg_D1 <- as.matrix(cbind( pars[10] + pars[11]*tri_var$scalex1 + pars[12]*tri_var$scalez1+
                              pars[13]*tri_var$z2, rep(pars[14],n), rep(pars[15],n)))
marginals1 <- list(marg_CP[1,], marg_WS[1,], marg_D1[1,]) 
marginals2 <- list(marg_CP[2,], marg_WS[2,], marg_D1[2,]) 


#Simulating the posterior predictive for the Group 1 storms
set.seed(99399)
for(i in 1:num_simulations){
vars_tri1[i,] <- rmvevd(1, dep = pars[16], model = "log", d = 3, mar = marginals1)
vars_tri2[i,] <- rmvevd(1, dep = pars[16], model = "log", d = 3, mar = marginals2)
}


quants <- c(0.025, 0.975)
apply(vars_tri1 , 2 , quantile , probs = quants , na.rm = TRUE )
apply(vars_tri2 , 2 , quantile , probs = quants , na.rm = TRUE )

true_vals1 <- c(tri_var$z1[1], tri_var$x1[1], tri_var$x2[1])
true_vals1

true_vals2 <- c(tri_var$z1[2], tri_var$x1[2], tri_var$x2[2])
true_vals2

## Finding the percentile for the true-values
ecdf_fun <- function(x,perc) ecdf(x)(perc)
alpha_minCP2016_1 <- ecdf_fun(vars_tri1[,1],true_vals1[1])
2*min(alpha_minCP2016_1, 1-alpha_minCP2016_1)

alpha_minCP2016_2 <- ecdf_fun(vars_tri2[,1],true_vals2[1])
2*min(alpha_minCP2016_2, 1-alpha_minCP2016_2)

alpha_maxWS2016_1 <- ecdf_fun(vars_tri1[,2],true_vals1[2])
2*min(alpha_maxWS2016_1, 1- alpha_maxWS2016_1)

alpha_maxWS2016_2 <- ecdf_fun(vars_tri2[,2],true_vals2[2])
2*min(alpha_maxWS2016_2, 1- alpha_maxWS2016_2)

alpha_damage2016_1 <- ecdf_fun(vars_tri1[,3], true_vals1[3])
2*min(alpha_damage2016_1, 1-alpha_damage2016_1)

alpha_damage2016_2 <- ecdf_fun(vars_tri2[,3], true_vals2[3])
2*min(alpha_damage2016_2, 1-alpha_damage2016_2)



#2016 Simulated minCP with observed value for the two hurricanes
pdf("Pred2016_minCP_moda_triEVD_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(vars_tri1[,1],freq=F,xlab='log min central Pressure',main='Hurricane Hermine 2016',ylim=c(0,0.7), axes = F, breaks=c(seq(-4,7.5,by=0.5)),cex.lab=1.25)
axis(side=2,at=seq(0,0.7,0.15),pos=-4,cex.axis=1.25)
axis(side=1,at=seq(-4,7.5,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-4,7.5,0.5),pos=0.7,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.7,0.02),pos=7.5,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$z1[1]+0.01,y0=0, x1= tri_var$z1[1] +0.01,y1= 0.52,col='red',lwd=3,lty=3)
segments(x0= tri_var$z1[1] -0.01,y0=0, x1= tri_var$z1[1] -0.01,y1= 0.52,col='red',lwd=3,lty=3)
segments(x0=-4,y0=0.60,x1=-4,y1=0.70)
dev.off()

pdf("Pred2016Matthew_minCP_moda_triEVD_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(vars_tri2[,1],freq=F,xlab='log min central Pressure',main='Hurricane Matthew 2016',ylim=c(0,0.7), axes = F, breaks=c(seq(-4,7.5,by=0.5)),cex.lab=1.25)
axis(side=2,at=seq(0,0.7,0.15),pos=-4,cex.axis=1.25)
axis(side=1,at=seq(-4,7.5,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-4,7.5,0.5),pos=0.7,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.7,0.02),pos=7.5,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$z1[2]+0.01,y0=0, x1= tri_var$z1[2]+0.01,y1= 0.58,col='red',lwd=3,lty=3)
segments(x0= tri_var$z1[2]-0.01,y0=0, x1= tri_var$z1[2]-0.01,y1= 0.58,col='red',lwd=3,lty=3)
segments(x0=-4,y0=0.60,x1=-4,y1=0.70)
dev.off()

#2016 Simulated WS for hurricance Hermine in 2016
pdf("Pred2016Hermine_WS_moda_triEVD_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(vars_tri1[,2],freq=F,xlab='log max WS',main='Hurricane Hermine 2016', xlim = c(0,5), ylim = c(0,1.5), axes = F)
axis(side=2,at=seq(0,1.5,0.2),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,5,1),pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,5,0.5),pos=1.5,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,1.5,0.01),pos=5,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$x1[1]+0.01,y0=0, x1= tri_var$x1[1]+0.01,y1= 1.1,col='red',lwd=3,lty=3)
segments(x0= tri_var$x1[1]-0.01 ,y0=0, x1=tri_var$x1[1]-0.01,y1= c(1.1, 1.1),col='red',lwd=3,lty=3)
segments(x0=0,y0=1.4,x1=0,y1=1.5)
dev.off()


pdf("Pred2016Matthew_WS_moda_triEVD_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(vars_tri2[,2],freq=F,xlab='log max WS',main='Hurricane Hermine 2016', ylim = c(0,1.5), xlim = c(0,6), axes = F)
axis(side=2,at=seq(0,1.5,0.2),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,6,1),pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,6,0.5),pos=1.5,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,1.5,0.01),pos=6,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$x1[2]+0.01,y0=0, x1= tri_var$x1[2]+0.01,y1= 0.85,col='red',lwd=3,lty=3)
segments(x0= tri_var$x1[2]-0.01 ,y0=0, x1=tri_var$x1[2]-0.01,y1= 0.85,col='red',lwd=3,lty=3)
segments(x0=0,y0=1.4,x1=0,y1=1.5)

dev.off()


#2016 Simulated Damages with observed value
pdf("Pred2016Hermine_Damages_moda_triEVD_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(vars_tri1[,3],freq=F,xlab='log max WS',main='Hurricane Matthew 2016', xlim = c(0,25), ylim = c(0, 0.3), axes = F)
axis(side=2,at=seq(0,0.3,0.01),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,25,5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,25,0.5),pos=0.3,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.3,0.02),pos=25,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$x2[1] + 0.01,y0=0, x1= tri_var$x2[1] + 0.01,y1= c(0.13, 0.13),col='red',lwd=3,lty=3)
segments(x0= tri_var$x2[1] - 0.01,y0=0, x1= tri_var$x2[1] - 0.01,y1= c(0.13, 0.13),col='red',lwd=3,lty=3)
#segments(x0=0,y0=1.4,x1=0,y1=1.5)
dev.off()


#2016 Simulated Damages with observed value
pdf("Pred2016Matthew_Damages_moda_triEVD_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(vars_tri2[,3],freq=F,xlab='log max WS',main='Hurricane Matthew 2016', xlim = c(0,28), ylim = c(0, 0.3), axes = F)
axis(side=2,at=seq(0,0.3,0.01),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,28,5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,28,0.5),pos=0.3,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.3,0.02),pos=28,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$x2[2] + 0.01,y0=0, x1= tri_var$x2[2] + 0.01,y1= 0.115,col='red',lwd=3,lty=3)
segments(x0= tri_var$x2[2] - 0.01,y0=0, x1= tri_var$x2[2] - 0.01,y1= 0.115,col='red',lwd=3,lty=3)
#segments(x0=0,y0=1.4,x1=0,y1=1.5)
dev.off()

