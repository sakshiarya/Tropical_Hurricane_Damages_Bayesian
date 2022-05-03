##############################################################################
##############################################################################
#####  Title: Posterior Predictives for Tri-variate hierarchical model #######
##(with half-normal damages) using Best Estimate for 2016 based on 1960-2015##
#Author: Sakshi Arya (aryax010@umn.edu)
#Objective: Create models and run MCMC
#Created: 08/27/21
#Last updated: 08/27/21
##############################################################################
##############################################################################

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
library(tseries)
library(stringr)
library(pscl)
library(evd)
library(extRemes)

#Run all data creation and hyperparameter code
filep<-"/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/BayesianBivEVDMCMC/Predictions/Predicting 2016/HierWithHalfNormalDamages/"
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
tri_var$x2 <- tri_var$CURRENT.DAMAGE.2021
tri_var$z1 <- tri_var$logminCP
tri_var$scalez1 <- scale(tri_var$logminCP)
tri_var$z2 <- scale(tri_var$avgLat)

# #Annualized Covariates (May/June Averages)
# Annual_Cov1 <- read.csv(file.path(filedata,'Derived Data Sets 2021/Annual_Covariates_1960_2019.csv'))
# Annual_Cov <- Annual_Cov1[Annual_Cov1$Year>=2016,]

#testset<- as.numeric(Annual_Cov[Annual_Cov$Year == 2016,4:9])

#Loading in the MCMC Chains
load(paste(filep,'hierNS_HFNDamages_1e6_27August_1960to2015.RData',sep=''))



#Saving the hierarchical estimates (means of the posterior)
pars <- colMeans(results)


#Making vectors to store simulated posterior predictive
#lambdas <-vector()
num_simulations<-10^5
n <- dim(storm_damage_2016)[1]
WS <-matrix(NA, ncol = n, nrow = num_simulations)
CP <-matrix(NA, ncol = n, nrow = num_simulations)
D1 <-matrix(NA, ncol = n, nrow = num_simulations)



#Simulating the posterior predictive for the Group 1 storms
set.seed(99399)
for(i in 1:num_simulations){
  CP[i,] <- revd(n, loc = pars[1] + pars[2]*tri_var$z2, scale = rep(pars[3],n), 
               shape = rep(pars[4], n), type = "GEV")
  WS[i,] <- revd(n, loc = pars[5] + pars[6]*tri_var$scalez1 + pars[7]*tri_var$z2,
                 scale = rep(pars[8],n), 
                 shape = rep(pars[9], n), type = "GEV")
  D1[i,] <- rlnorm(n, meanlog = pars[10] + pars[11]*tri_var$scalex1 + pars[12]*tri_var$scalez1+
                     pars[13]*tri_var$z2, sdlog = pars[14])
}

quants <- c(0.025, 0.975)
apply(CP , 2 , quantile , probs = quants , na.rm = TRUE )
apply(WS, 2 , quantile , probs = quants , na.rm = TRUE )
apply(log(D1), 2 , quantile , probs = quants , na.rm = TRUE )


true_vals1 <- c(tri_var$z1[1], tri_var$x1[1], log(tri_var$x2)[1])
true_vals1
true_vals2 <- c(tri_var$z1[2], tri_var$x1[2], log(tri_var$x2)[2])
true_vals2

## Finding the percentile for the true-values
ecdf_fun <- function(x,perc) ecdf(x)(perc)
alpha_minCP2016_1 <- ecdf_fun(CP[,1],true_vals1[1])
2*min(alpha_minCP2016_1, 1-alpha_minCP2016_1)

alpha_minCP2016_2 <- ecdf_fun(CP[,2],true_vals2[1])
2*min(alpha_minCP2016_2, 1-alpha_minCP2016_2)

alpha_maxWS2016_1 <- ecdf_fun(WS[,1],true_vals1[2])
2*min(alpha_maxWS2016_1, 1- alpha_maxWS2016_1)

alpha_maxWS2016_2 <- ecdf_fun(WS[,2],true_vals2[2])
2*min(alpha_maxWS2016_2, 1- alpha_maxWS2016_2)

alpha_damage2016_1 <- ecdf_fun(log(D1[,1]), true_vals1[3])
2*min(alpha_damage2016_1, 1-alpha_damage2016_1)

alpha_damage2016_2 <- ecdf_fun(log(D1[,2]), true_vals2[3])
2*min(alpha_damage2016_2, 1-alpha_damage2016_2)

## Finding the percentile for the true-values
ecdf_fun <- function(x,perc) ecdf(x)(perc)
alpha_minCP2016_1 <- ecdf_fun(CP[,1],true_vals1[1])
2*min(alpha_minCP2016_1, 1-alpha_minCP2016_1)

alpha_minCP2016_2 <- ecdf_fun(CP[,2],true_vals2[1])
2*min(alpha_minCP2016_2, 1-alpha_minCP2016_2)

alpha_maxWS2016_1 <- ecdf_fun(WS[,1],true_vals1[2])
2*min(alpha_maxWS2016_1, 1- alpha_maxWS2016_1)

alpha_maxWS2016_2 <- ecdf_fun(WS[,2],true_vals2[2])
2*min(alpha_maxWS2016_2, 1- alpha_maxWS2016_2)

alpha_damage2016_1 <- ecdf_fun(D1[,1], true_vals1[3])
2*min(alpha_damage2016_1, 1-alpha_damage2016_1)

alpha_damage2016_2 <- ecdf_fun(D1[,2], true_vals2[3])
2*min(alpha_damage2016_2, 1-alpha_damage2016_2)

#2016 Simulated Frequency for TS-2 with observed value
pdf("PredMatthew2016_minCP_HNdamagehier_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(CP[,2],freq=F,xlab='log of Minimum Central Pressure',main='',ylim=c(0,0.7), axes = F, breaks=c(seq(-4,7.5,by=0.5)),cex.lab=1.25)
axis(side=2,at=seq(0,0.7,0.15),pos=-4,cex.axis=1.25)
axis(side=1,at=seq(-4,7.5,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-4,7.5,0.5),pos=0.7,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.7,0.02),pos=7.5,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$z1[2] + 0.01,y0=0, x1= tri_var$z1[2] + 0.01,y1= 0.6,col='red',lwd=3,lty=3)
segments(x0= tri_var$z1[2] - 0.01,y0=0, x1= tri_var$z1[2] - 0.01,y1= 0.6,col='red',lwd=3,lty=3)
segments(x0=-4,y0=0.60,x1=-4,y1=0.70)
dev.off()

#2016 Simulated Landfall for TS-2 with observed value
pdf("PredMatthew2016_WS_HNdamagehier_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(WS[,2],freq=F,xlab='log of Maximum Wind Speed',main='',ylim=c(0,1.2), axes = F, breaks=c(seq(0,7,by=0.25)),cex.lab=1.25)
axis(side=2,at=seq(0,1.2,0.1),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,7,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,7,0.2),pos=1.2,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,1.2,0.02),pos=7,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= tri_var$x1[2] + 0.01,y0=0, x1= tri_var$x1[2] + 0.01,y1= 1.02 ,col='red',lwd=3,lty=3)
segments(x0= tri_var$x1[2] - 0.01,y0=0, x1= tri_var$x1[2] - 0.01,y1= 1.02,col='red',lwd=3,lty=3)

dev.off()

#2011 Simulated Damages for TS-2 with observed value
pdf("Pred2016_Damages_methodc_HN_27Aug.pdf")
par(mar=c(3,3,0,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(log(D1[,2]),freq=F,xlab='log of Storm Damage',main='',ylim=c(0,0.27), axes = F, breaks=c(0,1,2:35),cex.lab=1.25)
axis(side=2,at=seq(0,0.27,0.05),pos=0,cex.axis=1.25)
axis(side=1,at=seq(0,35,5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(0,35,0.2),pos=0.27,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.27,0.02),pos=35,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= log(storm_damage_2016$CURRENT.DAMAGE.2021 +0.1),y0=0, x1= log(storm_damage_2016$CURRENT.DAMAGE.2021 +0.1),y1= c(0.19, 0.10),col='red',lwd=3,lty=3)
segments(x0= log(storm_damage_2016$CURRENT.DAMAGE.2021 -0.1),y0=0, x1= log(storm_damage_2016$CURRENT.DAMAGE.2021 -0.1),y1= c(0.19, 0.10),col='red',lwd=3,lty=3)
segments(x0=0,y0=0.20,x1=0,y1=0.22)
dev.off()





################## Prediction method 2 #################
dim(results)
results_sub <- results[(9*(10^5)+1):10^6,]
num_simulations <- 100
post_samples <- dim(results_sub)[1]

WS <-matrix(NA, ncol = post_samples, nrow = num_simulations)
CP <-matrix(NA, ncol = post_samples, nrow = num_simulations)
D1 <-matrix(NA, ncol = post_samples, nrow = num_simulations)

start_time <- Sys.time()
for(j in 1:post_samples){
  CP[,j] <- revd(num_simulations, loc = pars[1] + pars[2]*tri_var$z2[1], scale = pars[3], 
                 shape = pars[4], type = "GEV")
  WS[,j] <- revd(num_simulations, loc = pars[5] + pars[6]*tri_var$z1[1] + pars[7]*tri_var$z2[1],
                 scale = pars[8], 
                 shape = pars[9], type = "GEV")
  D1[,j] <- rlnorm(num_simulations, meanlog = pars[10] + pars[11]*tri_var$scalex1[1] + pars[12]*tri_var$z1[1]+
                     pars[13]*tri_var$z2[1], sdlog = pars[14])
  print(j)
}


## Plot for damages
values<-data.frame(val=density(log(D1[,3]),bw=0.25)$x)
pdf_func <- function(i){
  return(density(log(D1[, i]),bw=0.25)$y)
}
result_D1 <- lapply(1:post_samples, pdf_func)
#result_D1 <- lapply(1:poster, pdf_func)
result_Means_D1 <- rowMeans(do.call(cbind, result_D1), na.rm = TRUE)
result_Mean_D1 <- cbind(values, result_Means_D1)
head(result_Mean_D1)

pdf("Pred2016_Damages_methodc_predb_HN_27Aug.pdf")
par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(cbind(result_Mean_D1[,1],result_Mean_D1[,2]),type='h',xlab='TS-2 Storm Damage',main='',ylab='Density', cex.lab = 1.25, cex.axis = 1.25, lwd = 1.5)
damage.truth.2016 <- log(testsettruth[testsettruth$Year==2016 & testsettruth$Cat_HURDAT=='TS-2',]$damage)
diff.ts2.damage.from.truth.2016 <- abs(ts2.damage.truth.2016 - result_Mean_D1[,1])
segments(x0= ts2.damage.truth.2016-0.1,y0=0, x1= ts2.damage.truth.2016-0.1,y1= result_Mean_D1[which(diff.ts2.damage.from.truth.2016 == min(diff.ts2.damage.from.truth.2016)),2],col='red',lwd=3,lty=2)
segments(x0= ts2.damage.truth.2016+0.1,y0=0, x1= ts2.damage.truth.2016+0.1,y1= result_Mean_D1[which(diff.ts2.damage.from.truth.2016 == min(diff.ts2.damage.from.truth.2016)),2],col='red',lwd=3,lty=2)


## Plot for minCP
values<-data.frame(val=density(CP[,3],bw=0.25)$x)
pdf_func <- function(i){
  return(density(CP[, i],bw=0.25)$y)
}
result_CP <- lapply(1:post_samples, pdf_func)
#result_CP <- lapply(1:poster, pdf_func)
result_Means_CP <- rowMeans(do.call(cbind, result_CP), na.rm = TRUE)
result_Mean_CP <- cbind(values, result_Means_CP)
head(result_Mean_CP)

par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(cbind(result_Mean_CP[,1],result_Mean_CP[,2]),type='h',xlab='TS-2 Storm Damage',main='',ylab='Density', cex.lab = 1.25, cex.axis = 1.25, lwd = 1.5, xlim = c(-3,7))
#text(5,0.15,paste(round(avg.zeros,2),'chance of\n $0 Damage'), cex = 1.15)
#text(25.4,0.17,paste(round(1-avg.zeros,2),'chance of\n log(Damage)\n in this distribution'), cex = 1.12)
damage.truth.2016 <- log(storm_damage_2016$minCP)
diff.damage.from.truth.2016 <- abs(damage.truth.2016 - result_Mean_CP[,1])
segments(x0= damage.truth.2016-0.1,y0=0, x1= damage.truth.2016-0.1,y1= c(0.4, 0.4),col='red',lwd=3,lty=2)
segments(x0= damage.truth.2016+0.1,y0=0, x1= damage.truth.2016+0.1,y1= c(0.4,0.4),col='red',lwd=3,lty=2)



### Wind Speed plots
values<-data.frame(val=density(WS[,3],bw=0.25)$x)
pdf_func <- function(i){
  return(density(WS[, i],bw=0.25)$y)
}
result_WS <- lapply(1:post_samples, pdf_func)
#result_WS <- lapply(1:poster, pdf_func)
result_Means_WS <- rowMeans(do.call(cbind, result_WS), na.rm = TRUE)
result_Mean_WS <- cbind(values, result_Means_WS)
head(result_Mean_WS)

par(mar=c(4,4,0,0)+0.1,mgp=c(2.5, 0.8, 0))#sets margins of plotting area
plot(cbind(result_Mean_WS[,1],result_Mean_WS[,2]),type='h',xlab='TS-2 Storm WS',main='',ylab='Density', cex.lab = 1.25, cex.axis = 1.25, lwd = 1.5)
WS.truth.2016 <- log(storm_damage_2016$maxWS)
diff.WS.from.truth.2016 <- cbind(abs(WS.truth.2016[1] - result_Mean_WS[,1]), abs(WS.truth.2016[2] - result_Mean_WS[,1]))

segments(x0= log(storm_damage_2016$maxWS -0.1),y0=0, x1= log(storm_damage_2016$maxWS -0.1),y1= result_Mean_WS[c(which.min(diff.WS.from.truth.2016[,1]),
                                                                                                                which.min(diff.WS.from.truth.2016[,2])),2],col='red',lwd=3,lty=2)
segments(x0= log(storm_damage_2016$maxWS + 0.1),y0=0, x1= log(storm_damage_2016$maxWS +0.1),y1= result_Mean_WS[c(which.min(diff.WS.from.truth.2016[,1]),
                                                                                                                 which.min(diff.WS.from.truth.2016[,2])),2],col='red',lwd=3,lty=2)







