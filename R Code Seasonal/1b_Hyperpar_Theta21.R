#################################################################
#################################################################
#Title: Hyperparameters for Landfall
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Determine hurricane landfall hyperparameters, theta
#Created: 2/25/16
#Last updated: 3/18/2021 (by Sakshi Arya)
#################################################################
#################################################################
#rm(list = ls())
#Loads all necessary libraries
library(lme4)
library(car)
library(ggplot2)

#################################################################
#################################################################
#Reading in Derived Data Sets
#################################################################
#################################################################

filedata<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Data Sets/'

#Storm-wise data
Storm_Data <-read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Storm_1960_2019.csv'))

#Annualized data
Annual_Data<- read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Annual_1960_2019.csv'))

#################################################################
#################################################################
# Examining whether there is Clustering of Landfall within a year
#################################################################
#################################################################
#Putting storm data in time order
clustering<-Storm_Data[order(Storm_Data$Year, Storm_Data$StartMonth),]
clustering$DamageInd <- as.factor(clustering$DamageInd)
clustering$DamageInd<-relevel(clustering$DamageInd,'No')
clustering$DI<-as.numeric(clustering$DamageInd)-1

clustering$lagDI<-NA

#Creating a 1 period lag for damage indicator
clustering[2:dim(clustering)[1],]$lagDI<-clustering[1:(dim(clustering)[1]-1),]$DI

#Assigning the first storm of the year no lag value
for(year in 1960:2019){
  clustering[clustering$Year==year,][1,]$lagDI<-NA
}

#Logit Normal Mixed Model with Fixed lag effect + Random int by year
summary(mod.re<-glmer(DI ~ lagDI +(1|Year),data= clustering,family='binomial', control=glmerControl(optimizer="bobyqa")))
Anova(mod.re)

#Fixed lag effect only
summary(mod.fe<-glm(DI ~ lagDI, data= clustering,family='binomial'))
Anova(mod.fe)

#LRT for variance component
anova(mod.re, mod.fe)

############################################################
############################################################
#Hyperparameters for Theta
# Mean of beta distribution:
# alpha/(alpha+beta) = Damage_Prob
# (Damage_Prob/(1-Damage_Prob)) alpha= beta

# Approx. Median of beta distribution:
# ~(alpha-1/3)/(alpha+beta-2/3) = Damage_Prob
#  ((1-Damage_Prob)* alpha -1/3)/Damage_Prob+2/3 = beta
#################################################################
#################################################################
#Changing factor level orders for graphing
Storm_Data$DamageInd <- as.factor(Storm_Data$DamageInd)
Storm_Data$Cat_HURDAT <- as.factor(Storm_Data$Cat_HURDAT)
Annual_Data$Cat_HURDAT <- as.factor(Annual_Data$Cat_HURDAT)
Storm_Data$DamageInd<-relevel(Storm_Data$DamageInd,'No')
Storm_Data$Cat_HURDAT <-relevel(Storm_Data$Cat_HURDAT,'TS-2')
Annual_Data$Cat_HURDAT <-relevel(Annual_Data$Cat_HURDAT,'TS-2')

#pdf("Proportion_Landfall_new.pdf")
ggplot(Storm_Data,aes(x= Cat_HURDAT,fill= DamageInd))+ geom_bar(position='fill',alpha=0.75)+xlab('Maximum Storm Category')+ylab('Proportion Landfall/No Landfall (1960-2013)')+ scale_fill_manual('Landfall Indicator',values=c("grey","blue"))+facet_wrap(~Year)+ theme_classic()
#dev.off()

#Individual storms
Damage_Prob_Storm <-with(Storm_Data, tapply(DamageInd=='Yes', Cat_HURDAT,mean))
Damage_Prob_Annual <-with(Annual_Data, tapply(Landfall/Freq, Cat_HURDAT,function(x){mean(x,na.rm=T)}))

#Settings for Beta Distribution Prior
alpha =c(1,1)

#Mean setting
#beta = alpha * Damage_Prob/(1-Damage_Prob)

#Median setting (makes more sense based on graph)
beta_storm =((1-Damage_Prob_Storm)* alpha -1/3)/Damage_Prob_Storm +2/3
beta_annual =((1-Damage_Prob_Annual)* alpha -1/3)/Damage_Prob_Annual +2/3

#Using annual 
hyper_Theta<-cbind(alpha, beta_annual)

#Plotting the output of the priors with the hyperparameters specified
#par(mar=c(3,3,0,0)+0.1, mgp =c(2,1,0))#sets margins of plotting area
#pdf("Prior_TS2_Landfall_Theta_new.pdf")
curve(dbeta(x, alpha[1], beta_annual[1]),xlab= expression (paste(theta[1])),ylab='Density',from=0,to=1,cex.axis=1.25, cex.lab=1.5)
#dev.off()

#pdf("Prior_35_Landfall_Theta_new.pdf")
curve(dbeta(x, alpha[2], beta_annual[2]),xlab= expression (paste(theta[2])),ylab='Density',from=0,to=1,cex.axis=1.25, cex.lab=1.5)
#dev.off()
