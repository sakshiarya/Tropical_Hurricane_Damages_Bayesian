#################################################################
#################################################################
#Title: Posterior Predictives for Tri-variate EVD model using Best Estimate for 2022 based on 1960-2021
#Author: Sakshi Arya (aryax010@umn.edu)
#Objective: Create models and run MCMC
#Created: 08/27/21
#Last updated: 03/24/2024
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

## File storage destinations
filep <- "/Users/ska5950/Library/CloudStorage/GoogleDrive-arya.sakshi44@gmail.com/My\ Drive/Research/2022/Ansu/BayesStorm/Tropical_Hurricane_Damages_Bayesian/2024_retake/R\ Code/Predictions"
filefigs <- "/Users/ska5950/Library/CloudStorage/GoogleDrive-arya.sakshi44@gmail.com/My\ Drive/Research/2022/Ansu/BayesStorm/Tropical_Hurricane_Damages_Bayesian/2024_retake/R\ Code/Predictions"
file.derived.22 <- "/Users/ska5950/Library/CloudStorage/GoogleDrive-arya.sakshi44@gmail.com/My\ Drive/Research/2022/Ansu/BayesStorm/Tropical_Hurricane_Damages_Bayesian/2024_retake/Data\ 2022"


## Reading the data:
storms<-read.csv(file.path(file.derived.22,'Categorized_Storm_1960_2022.csv'))
data_cov <- read.delim(file.path(file.derived.22,'Annual_Covariates_1960_2022.csv'), header = TRUE, sep = ",")
head(data_cov)
data_cov <- data_cov[,-which(colnames(data_cov) == "YR")]



names(storms)
storm_damage<-storms[storms$minCP< Inf & storms$CURRENT.DAMAGE.2022>0 & !is.na(storms$CURRENT.DAMAGE.2022),]
head(storm_damage)
tail(storm_damage)

data_combined <- left_join(storm_damage, data_cov, by = "Year")
head(data_combined)
dim(data_combined)


tri_var<-data_combined[,c('minCP','maxWS','CURRENT.DAMAGE.2022')]
tri_var$logdamage<-log(tri_var$CURRENT.DAMAGE.2022)
tri_var$logmaxWS<-log(tri_var$maxWS)

tri_var$logminCP<-log(1013-tri_var$minCP)
tri_var$avgLat <- storm_damage$avgLat
tri_var$avgLon <- storm_damage$avgLon
tri_var$StartMonth <- storm_damage$StartMonth
tri_var$Year <- storm_damage$Year
tri_var$max_Cat <- storm_damage$max_Cat
dim(tri_var)
tri_var$x1 <- tri_var$logmaxWS
tri_var$x2 <- tri_var$logdamage
tri_var$z1 <- tri_var$logminCP
tri_var[15:21] <- data_combined[,17:23]

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
tri_var <- tri_var[,-which(colnames(tri_var) %in% c("max_Cat", "total"))]

tri_var_sub1 <- tri_var[1:10]

tri_var_year_comb <- data.frame(cbind(tri_var, year= data_combined$Year))
covs_minCP <- colnames(tri_var_sub1[,-which(colnames(tri_var_sub1) %in% c("NAO","SOI", "ANOM.3.4","Sunspots"))])
covs_maxWS <- "scalez1"
covs_damages <- c("scalez1", "scalex1",colnames(tri_var_sub1[,-which(colnames(tri_var_sub1) %in% c("avgLon","NAO","Sunspots"))]))
testsettruth<- tri_var_year_comb[tri_var_year_comb$year==2022,-which(colnames(tri_var_year_comb) %in% "year")]
#densityPlot(tri_var2015$logdamage)


# #Annualized Covariates (May/June Averages)
# Annual_Cov1 <- read.csv(file.path(filedata,'Derived Data Sets 2021/Annual_Covariates_1960_2019.csv'))
# Annual_Cov <- Annual_Cov1[Annual_Cov1$Year>=2021,]

#testset<- as.numeric(Annual_Cov[Annual_Cov$Year == 2021,4:9])

#Loading in the MCMC Chains
load(file.path(filep,"hierNSfulldat_1e6_1960to2021.RData"))



#Saving the hierarchical estimates (means of the posterior)
pars <- colMeans(results)


#Making vectors to store simulated posterior predictive
#lambdas <-vector()
num_simulations<-10^5
n <- dim(testsettruth)[1]
d <- (length(covs_minCP) + length(covs_maxWS) + length(covs_damages) + 3)
WS <-matrix(NA, ncol = n, nrow = num_simulations)
CP <-matrix(NA, ncol = n, nrow = num_simulations)
D1 <-matrix(NA, ncol = n, nrow = num_simulations)


set.seed(99399)
for(i in 1:num_simulations){
  CP[i,] <- revd(n, loc = (pars[1] + colSums(t(testsettruth[,covs_minCP])*(pars[2:(length(covs_minCP)+1)]))), scale = rep(pars[d+2],n), 
                 shape =  rep(pars[d + 1], n), type = "GEV")
  WS[i,] <- revd(n, loc = pars[length(covs_minCP) + 2] + colSums((t(testsettruth[, covs_maxWS])*pars[length(covs_minCP) + length(covs_maxWS) + 2])),
                 scale = rep(pars[d+4],n), 
                 shape = rep(pars[d+3],n), type = "GEV")
  D1[i,] <- revd(n, loc = pars[length(covs_minCP) + length(covs_maxWS) + 3] + 
                   colSums(t(testsettruth[,covs_damages]*pars[(length(covs_minCP) + length(covs_maxWS) + 4):(length(covs_minCP) + length(covs_maxWS) + length(covs_damages) + 3)])),
                 scale = rep(pars[d+6],n), 
                 shape = rep(pars[d+5],n), type = "GEV")
}

quants <- c(0.025, 0.975)
apply(CP , 2 , quantile , probs = quants , na.rm = TRUE )
apply(WS, 2 , quantile , probs = quants , na.rm = TRUE )
apply(D1, 2 , quantile , probs = quants , na.rm = TRUE )


true_vals1 <- c(testsettruth$z1[1], testsettruth$x1[1], testsettruth$x2[1])
true_vals1
true_vals2 <- c(testsettruth$z1[2], testsettruth$x1[2], testsettruth$x2[2])
true_vals2
true_vals3 <- c(testsettruth$z1[3], testsettruth$x1[3], testsettruth$x2[3])
true_vals3

## Finding the percentile for the true-values
ecdf_fun <- function(x,perc) ecdf(x)(perc)
alpha_minCP2021_1 <- ecdf_fun(CP[,1],true_vals1[1])
2*min(alpha_minCP2021_1, 1-alpha_minCP2021_1)

alpha_minCP2021_2 <- ecdf_fun(CP[,2],true_vals2[1])
2*min(alpha_minCP2021_2, 1-alpha_minCP2021_2)

alpha_minCP2021_3 <- ecdf_fun(CP[,3],true_vals3[1])
2*min(alpha_minCP2021_3, 1-alpha_minCP2021_3)

alpha_maxWS2021_1 <- ecdf_fun(WS[,1],true_vals1[2])
2*min(alpha_maxWS2021_1, 1- alpha_maxWS2021_1)

alpha_maxWS2021_2 <- ecdf_fun(WS[,2],true_vals2[2])
2*min(alpha_maxWS2021_2, 1- alpha_maxWS2021_2)

alpha_damage2021_1 <- ecdf_fun(D1[,1], true_vals1[3])
2*min(alpha_damage2021_1, 1-alpha_damage2021_1)

alpha_damage2021_2 <- ecdf_fun(D1[,2], true_vals2[3])
2*min(alpha_damage2021_2, 1-alpha_damage2021_2)



#2021 Simulated Frequency for TS-2 with observed value
pdf(file.path(filefigs,"Pred2017Harvey_minCP_methodb_hierEVD_27April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist((CP[,1]),freq=F,xlab='log of Minimum Central Pressure',main='Hurricane Harvey 2017',ylim=c(0,0.6), breaks=c(seq(-9,6,by=0.5)), axes = F,cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.6,0.15),pos=-9,cex.axis=1.25)
axis(side=1,at=seq(-9,6,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-9,6,0.5),pos=0.6,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.6,0.02),pos=6,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= (testsettruth$z1[1] +0.05),y0=0, x1= (testsettruth$z1[1] +0.05),y1= c(0.525, 0.525),col='red',lwd=3,lty=3)
segments(x0= testsettruth$z1[1] -0.05,y0=0, x1= testsettruth$z1[1] - 0.05,y1= c(0.525, 0.525),col='red',lwd=3,lty=3)
segments(x0=-9,y0=0.60,x1=-9,y1=0.60)
dev.off()

pdf(file.path(filefigs,"Pred2017Irma_minCP_modb_hierEVD_27April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(CP[,2],freq=F,xlab='log min central Pressure',main='Hurricane Irma 2017',ylim=c(0,0.65), axes = F, breaks=c(seq(-11,6,by=0.5)),cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.65,0.15),pos=-11,cex.axis=1.25)
axis(side=1,at=seq(-11,6,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-11,6,0.5),pos=0.65,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.65,0.02),pos=6,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$z1[2]+0.05,y0=0, x1= testsettruth$z1[2]+0.05,y1= 0.57,col='red',lwd=3,lty=3)
segments(x0= testsettruth$z1[2]-0.05,y0=0, x1= testsettruth$z1[2]-0.05,y1= 0.57,col='red',lwd=3,lty=3)
segments(x0=-11,y0=0.60,x1=-11,y1=0.65)
dev.off()

pdf(file.path(filefigs,"Pred2017Nate_minCP_modb_hierEVD_23April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(CP[,3],freq=F,xlab='log min central Pressure',main='Hurricane Nate 2017',ylim=c(0,0.6), axes = F, breaks=c(seq(-8,5.5,by=0.5)),cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.6,0.15),pos=-8,cex.axis=1.25)
axis(side=1,at=seq(-8,5.5,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-8,5.5,0.5),pos=0.6,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.6,0.02),pos=5.5,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$z1[3]+0.05,y0=0, x1= testsettruth$z1[3]+0.05,y1= 0.37,col='red',lwd=3,lty=3)
segments(x0= testsettruth$z1[3]-0.05,y0=0, x1= testsettruth$z1[3]-0.05,y1= 0.37,col='red',lwd=3,lty=3)
# segments(x0=-8,y0=0.60,x1=-8,y1=0.70)
dev.off()

#2021 Simulated WS for hurricance Hermine in 2021
pdf(file.path(filefigs,"Pred2017Harvey_WS_modb_hierEVD_23April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(WS[,1],freq=F,xlab='log max WS',main='Hurricane Harvey 2017', breaks=c(seq(-2,7,by=0.5)), ylim = c(0,0.5), axes = F,cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.5,0.2),pos=-2,cex.axis=1.25)
axis(side=1,at=seq(-2,7,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-2,7,0.5),pos=0.5,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.5,0.01),pos=7,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$x1[1]+0.05,y0=0, x1= testsettruth$x1[1]+0.05,y1= 0.435,col='red',lwd=3,lty=3)
segments(x0= testsettruth$x1[1]-0.05 ,y0=0, x1=testsettruth$x1[1]-0.05,y1= 0.435,col='red',lwd=3,lty=3)
segments(x0=-2,y0=0.4,x1=-2,y1=0.5)
dev.off()


pdf(file.path(filefigs,"Pred2017Irma_WS_modb_hierEVD_27April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(WS[,2],freq=F,xlab='log max WS',main='Hurricane Irma 2017', ylim = c(0,0.5), breaks = c(seq(-2,7,by=0.5)), axes = F, cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.5,0.05),pos=-2,cex.axis=1.25)
axis(side=1,at=seq(-2,7,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-2,7,0.5),pos=0.5,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.5,0.01),pos=7,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$x1[2]+0.05,y0=0, x1= testsettruth$x1[2]+0.05,y1= 0.435,col='red',lwd=3,lty=3)
segments(x0= testsettruth$x1[2]-0.05 ,y0=0, x1=testsettruth$x1[2]-0.05,y1= 0.435,col='red',lwd=3,lty=3)
# segments(x0=0,y0=1.4,x1=0,y1=1.5)

dev.off()


pdf(file.path(filefigs,"Pred2017Nate_WS_modb_hierEVD_27April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(WS[,3],freq=F,xlab='log max WS',main='Hurricane Nate 2017', ylim = c(0,0.5), breaks = c(seq(-2,7,by=0.5)), axes = F, cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.5,0.05),pos=-2,cex.axis=1.25)
axis(side=1,at=seq(-2,7,0.5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(-2,7,0.5),pos=0.5,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.5,0.01),pos=7,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$x1[3]+0.05,y0=0, x1= testsettruth$x1[3]+0.05,y1= 0.43,col='red',lwd=3,lty=3)
segments(x0= testsettruth$x1[3]-0.05 ,y0=0, x1=testsettruth$x1[3]-0.05,y1= 0.43,col='red',lwd=3,lty=3)
# segments(x0=0,y0=1.4,x1=0,y1=1.5)

dev.off()


#2017 Simulated Damages with observed value
pdf(file.path(filefigs,"Pred2017Harvey_Damages_moda_hierEVD_27April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(D1[,1],freq=F,xlab='log damages',main='Hurricane Harvey 2017', xlim = c(10,30), ylim = c(0, 0.25), axes = F, cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.25,0.05),pos=10,cex.axis=1.25)
axis(side=1,at=seq(10,30,5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(10,30,0.5),pos=0.25,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.25,0.05),pos=30,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$x2[1] + 0.05,y0=0, x1= testsettruth$x2[1] + 0.05,y1= c(0.03, 0.03),col='red',lwd=3,lty=3)
segments(x0= testsettruth$x2[1] - 0.05,y0=0, x1= testsettruth$x2[1] - 0.05,y1= c(0.03, 0.03),col='red',lwd=3,lty=3)
#segments(x0=0,y0=1.4,x1=0,y1=1.5)
dev.off()


#2017 Simulated Damages with observed value
pdf(file.path(filefigs,"Pred2017Irma_Damages_modb_hierEVD_27April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(D1[,2],freq=F,xlab='log damages',main='Hurricane Irma 2017', xlim = c(10,30), ylim = c(0, 0.25), axes = F, cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.25,0.05),pos=10,cex.axis=1.25)
axis(side=1,at=seq(10,30,5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(10,30,0.5),pos=0.25,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.25,0.05),pos=30,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$x2[2] + 0.05,y0=0, x1= testsettruth$x2[2] + 0.05,y1= 0.105,col='red',lwd=3,lty=3)
segments(x0= testsettruth$x2[2] - 0.05,y0=0, x1= testsettruth$x2[2] - 0.05,y1= 0.105,col='red',lwd=3,lty=3)
#segments(x0=0,y0=1.4,x1=0,y1=1.5)
dev.off()

#2017 Simulated Damages with observed value
pdf(file.path(filefigs,"Pred2017Nate_Damages_modb_hierEVD_27April23.pdf"))
par(mar=c(3,3,2,0),mgp=c(1.8, 0.8, 0))#sets margins of plotting area
hist(D1[,3],freq=F,xlab='log damages',main='Hurricane Nate 2017', xlim = c(10,30), ylim = c(0, 0.25), axes = F, cex.lab=1.5, cex.main = 1.5)
axis(side=2,at=seq(0,0.25,0.05),pos=10,cex.axis=1.25)
axis(side=1,at=seq(10,30,5),pos=0,cex.axis=1.25)
axis(side=3,at=seq(10,30,0.5),pos=0.25,lwd.ticks=0,labels=F)
axis(side=4,at=seq(0,0.25,0.05),pos=30,lwd.ticks=0,labels=F)
#text(6,0.16,paste(round(mean(D1==0),2),'chance of\n $0 Damage'), cex = 1.2)
#text(25,0.18,paste(round(1-mean(D1==0),2),'chance of\n log(Damage)\n in this distribution'), cex = 1.2)
segments(x0= testsettruth$x2[3] + 0.05,y0=0, x1= testsettruth$x2[3] + 0.05,y1= 0.205,col='red',lwd=3,lty=3)
segments(x0= testsettruth$x2[3] - 0.05,y0=0, x1= testsettruth$x2[3] - 0.05,y1= 0.205,col='red',lwd=3,lty=3)
#segments(x0=0,y0=1.4,x1=0,y1=1.5)
dev.off()








