#################################################################
#################################################################
#Title: Hyperparameters for Frequency
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Determine hurricane occurrence hyperparameters, beta
#			and phi
#Last Updated: 2/24/16
#Last edited: 3/18/2021 (by Sakshi (aryax010@umn.edu))
#################################################################
#################################################################
#rm(list = ls())
#Loads all necessary libraries
library(lme4)
# install.packages("DataCombine")
library(DataCombine)
library(car)
library(ggplot2)
library(MASS)
#################################################################
#################################################################
#Reading in Derived Data Sets
#################################################################
#################################################################

filedata<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/Data Sets'

#Storm-wise data
Storm_Data <-read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Storm_1960_2019.csv'))

#Annualized data
Annual_Data<- read.csv(file.path(filedata,'Derived Data Sets 2021/Categorized_Annual_1960_2019.csv'))

#Annualized Covariates (May/June Averages)
Annual_Cov <- read.csv(file.path(filedata,'Derived Data Sets 2021/Annual_Covariates_1960_2019.csv'))

#Years 1960-2019
GLM_Data <-merge(Annual_Data, Annual_Cov,by='Year')


################################################################################
#Category TS-2 Models
################################################################################
#Model 1 : Full Poisson TS-2
Poi1_all <-glm(Freq ~ NAO + SOI + AMO + ANOM.3.4 + Atl_SST + Sunspots, data= GLM_Data, 
               subset= Cat_HURDAT=='TS-2',family= poisson)
summary(Poi1_all)

f <- ~NAO + SOI + AMO + ANOM.3.4 + Atl_SST + Sunspots

#Model 2 : Intercept only Poisson TS-2
Poi1_int <-glm(Freq ~ 1,data= GLM_Data,
               subset= Cat_HURDAT=='TS-2', family= poisson)
summary(Poi1_int)

step(Poi1_int, scope = f, direction = "forward")
step(Poi1_all, scope = c(lower = ~ 1), direction = "backward" )

## Model 3: Selected Model Poisson
Poi1_select <- glm(Freq ~ NAO + AMO + ANOM.3.4 + Atl_SST, data = GLM_Data,subset= Cat_HURDAT=='TS-2', family= poisson )


#Model 4 : Full Neg. Bin. TS-2
NB1_all<-glm.nb(Freq ~ NAO + SOI + AMO + ANOM.3.4 + Atl_SST + Sunspots, data= GLM_Data,
                subset= Cat_HURDAT=='TS-2')
Phi1_all <-summary(NB1_all)$theta

xtable(summary(NB1_all)$coefficients[,1:2], digits = 3)

#Model 5 : Intercept only Neg. Bin. TS-2
NB1_int<-glm.nb(Freq ~ 1,data= GLM_Data,
                subset= Cat_HURDAT=='TS-2')
Phi1_int <-summary(NB1_int)$theta
summary(NB1_int)

## Stepwise selection
step(NB1_int, scope = f, direction = "forward")
step(NB1_all, scope = c(lower = ~ 1), direction = "backward" ) # Go with this one based on AIC value

## Model 6:
NB1_select <-glm.nb(Freq ~ NAO + AMO + ANOM.3.4 + Atl_SST, data= GLM_Data,
                subset= Cat_HURDAT=='TS-2')

#LRT Model 1 vs. Model 2
anova(Poi1_all, Poi1_int,test='Chisq')
# LRT Model 1 vs Model 3
anova(Poi1_all, Poi1_select,test='Chisq')
# LRT Model 2 vs Model 3
anova(Poi1_int, Poi1_select,test='Chisq')
# LRT Model 4 vs Model 5
anova(NB1_all, NB1_int, test = 'Chisq')
# LRT Model 4 vs Model 6
anova(NB1_all, NB1_select, test = 'Chisq')
#LRT Model 5 vs. Model 6
anova(NB1_select , NB1_int,test='Chisq')

#LRT Model 1 vs. Model 4
pchisq(2 * (logLik(NB1_all)[1] - logLik(Poi1_all)[1]), 
       df = abs(NB1_all$df.residual - Poi1_all$df.residual), lower.tail=FALSE)

#LRT Model 1 vs. Model 5
pchisq(2 * (logLik(Poi1_all)[1]-logLik(NB1_int)[1]), df = abs(Poi1_all$df.residual - NB1_int$df.residual), lower.tail=FALSE)

#LRT Model 1 vs. Model 6
pchisq(2 * (logLik(Poi1_all)[1]-logLik(NB1_select)[1]), df = abs(Poi1_all$df.residual - NB1_select$df.residual), lower.tail=FALSE)

# LRT Model 2 vs Model 4
pchisq(2 * (logLik(NB1_all)[1]-logLik(Poi1_all)[1]), df = abs(NB1_all$df.residual - Poi1_all$df.residual), lower.tail=FALSE)

# LRT Model 2 vs Model 5
pchisq(2 * (logLik(NB1_int)[1]-logLik(Poi1_int)[1]), df = abs(Poi1_int$df.residual - NB1_int$df.residual), lower.tail=FALSE)

# LRT Model 6 vs Model 2
pchisq(2 * (logLik(NB1_select)[1]-logLik(Poi1_int)[1]), df = abs(Poi1_int$df.residual - NB1_select$df.residual), lower.tail=FALSE)

#LRT Model 3 vs. Model 4
pchisq(2 * (logLik(NB1_all)[1] - logLik(Poi1_select)[1]), 
       df = abs(NB1_all$df.residual - Poi1_select$df.residual), lower.tail=FALSE)

#LRT Model 3 vs. Model 5
pchisq(2 * (logLik(Poi1_select)[1]- logLik(NB1_int)[1]), 
       df = abs(NB1_int$df.residual - Poi1_select$df.residual), lower.tail=FALSE)

#LRT Model 3 vs. Model 6
pchisq(2 * (logLik(NB1_select)[1] - logLik(Poi1_select)[1]), 
       df = abs(Poi1_select$df.residual - NB1_select$df.residual), lower.tail=FALSE)



#Lack of fit testing
pchisq(sum(resid(Poi1_select,type='pearson')^2), Poi1_select$df.residual, lower.tail=FALSE)
pchisq(sum(resid(Poi1_int,type='pearson')^2), Poi1_int$df.residual, lower.tail=FALSE)
pchisq(sum(resid(Poi1_select,type='pearson')^2), Poi1_select$df.residual, lower.tail=FALSE)
pchisq(sum(resid(NB1_all,type='pearson')^2),  NB1_all$df.residual, lower.tail=FALSE)
pchisq(sum(resid(NB1_int,type='pearson')^2),  NB1_int$df.residual, lower.tail=FALSE)
pchisq(sum(resid(NB1_select,type='pearson')^2),  NB1_select$df.residual, lower.tail=FALSE)

#Information Criteria for All TS-2 models
AIC(Poi1_all);BIC(Poi1_all)
AIC(Poi1_select); BIC(Poi1_select)
AIC(Poi1_int);BIC(Poi1_int)
AIC(NB1_all);BIC(NB1_all)
AIC(NB1_select); BIC(NB1_select)
AIC(NB1_int);BIC(NB1_int)

set.seed(1212)
test.data1<-data.frame(Value=c(rpois(1000, exp(coef(Poi1_int))),rnegbin(1000, exp(coef(NB1_int)), Phi1_int),GLM_Data[GLM_Data$Cat_HURDAT=='TS-2',]$Freq,rep(NA, 1000-length(GLM_Data[GLM_Data$Cat_HURDAT=='TS-2',]$Freq))), Type=c(rep('Pois. Sim.',1000),rep('NB Sim.',1000),rep('Empirical',1000)),check.names=F)

#pdf("DensityTS2_new.pdf")
ggplot(data=test.data1,aes(x= Value,y = ..density..,group=Type,fill=Type))+ stat_bin(position= position_dodge(4), binwidth = 5)+theme_bw()+ylab('Density')+scale_x_continuous(breaks = seq(2.5,39, by = 5),labels=c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39'),'# of Tropical Storm to Category 2 Tropical Cyclones')+ scale_fill_manual(values=c('grey0','grey50','grey75'))+theme(legend.position = c(0.7, 0.9), legend.justification = c(0, 1), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
#dev.off()
#ggplot(data=test.data1,aes(x= Value,y = ..density..,group=Type,fill=Type))+ stat_bin(position= position_dodge(4), binwidth = 5)+ labs(title = "Empirical Density of Tropical Storms - Cat. 2 Tropical Cyclones\n  vs. Poisson and NB Densities")+theme_bw()+ylab('Density')+scale_x_continuous(breaks = seq(2.5,39, by = 5),labels=c('0-4','5-9','10-14','15-19','20-24','25-29','30-34','35-39'),'# of Tropical Storm to Category 2 Tropical Cyclones')+ scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+theme(legend.position = c(0.7, 0.9), legend.justification = c(0, 1), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title=element_text(size=16),axis.text=element_text(size=14), legend.title=element_blank() ,legend.text=element_text(size=14), title =element_text(size=16))
################################################################################
#Category 3-5 Models
################################################################################

#Model 1 : Full Poisson 3-5
Poi2_all <-glm(Freq ~ NAO + SOI + AMO + ANOM.3.4 + Atl_SST + Sunspots,data= GLM_Data,subset=Cat_HURDAT=='3-5', family=poisson)
xtable(summary(Poi2_all)$coefficients[,1:2], digits = 3)

#Model 2 : Intercept + AMO Poisson 3-5
Poi2_int <-glm(Freq ~ 1,data= GLM_Data,subset=Cat_HURDAT=='3-5', family=poisson)

step(Poi2_int, scope = f, direction = "forward")
step(Poi2_all, scope = c(lower = ~ 1), direction = "backward" )

# Model 3 (selected): Adding Atl SST anyway because of practical reasons
Poi2_select <- glm(Freq ~ AMO + ANOM.3.4 + Atl_SST, data= GLM_Data,subset=Cat_HURDAT=='3-5', family=poisson)

#Model 4: Full NB 3-5
NB2_all <-glm.nb(Freq ~ NAO + SOI + AMO + ANOM.3.4 + Atl_SST + Sunspots, data= GLM_Data, subset=Cat_HURDAT=='3-5') # Not converging
#Model 5 : Intercept + AMO NB 3-5
NB2_int <-glm.nb(Freq ~ 1, data= GLM_Data, subset=Cat_HURDAT=='3-5')

step(NB2_int, scope = f, direction = "forward")


#Model 6 : Intercept + AMO NB 3-5
NB2_select <-glm.nb(Freq ~ AMO + ANOM.3.4,data= GLM_Data,subset=Cat_HURDAT=='3-5')



#LRT Model 1 vs. Model 2
anova(Poi2_all, Poi2_int,test='Chisq')
#LRT Model 1 vs. Model 3
anova(Poi2_all, Poi2_select,test='Chisq')
# LRT Model 2 vs. Model 3
anova(Poi2_int, Poi2_select,test='Chisq')

# LRT Model 4 vs. Model 5
anova(NB2_all, NB2_int, test = 'Chisq')
# LRT Model 4 vs. Model 6
anova(NB2_all, NB2_select, test = 'Chisq')
# LRT Model 5 vs. Model 6
anova(NB2_int, NB2_select, test = 'Chisq')

#LRT Model 1 vs. Model 4
pchisq(2 * (logLik(Poi2_all)[1]- logLik(NB2_all)[1]), 
       df = abs(NB2_all$df.residual - Poi2_all$df.residual), lower.tail=FALSE)

#LRT Model 1 vs. Model 5
pchisq(2 * (logLik(Poi2_all)[1]- logLik(NB2_int)[1]), 
       df = abs(NB2_int$df.residual - Poi2_all$df.residual), lower.tail=FALSE)

#LRT Model 1 vs. Model 6
pchisq(2 * (logLik(Poi2_all)[1]- logLik(NB2_select)[1]), 
       df = abs(NB2_select$df.residual - Poi2_all$df.residual), lower.tail=FALSE)


#LRT Model 2 vs. Model 4
pchisq(2 * (logLik(NB2_all)[1]- logLik(Poi2_int)[1]), 
       df = abs(NB2_all$df.residual - Poi2_int$df.residual), lower.tail=FALSE)

#LRT Model 2 vs. Model 5
pchisq(2 * (logLik(NB2_int)[1]- logLik(Poi2_int)[1]), 
       df = abs(NB2_int$df.residual - Poi2_int$df.residual), lower.tail=FALSE)

#LRT Model 2 vs. Model 6
pchisq(2 * (logLik(NB2_select)[1]- logLik(Poi2_int)[1]), 
       df = abs(NB2_select$df.residual - Poi2_int$df.residual), lower.tail=FALSE)

#LRT Model 3 vs. Model 4
pchisq(2 * (logLik(NB2_all)[1]- logLik(Poi2_select)[1]), 
       df = abs(NB2_all$df.residual - Poi2_select$df.residual), lower.tail=FALSE)

#LRT Model 3 vs. Model 5
pchisq(2 * (logLik(Poi2_select)[1]-logLik(NB2_int)[1]), 
       df = abs(NB2_int$df.residual - Poi2_select$df.residual), lower.tail=FALSE)

#LRT Model 3 vs. Model 6
pchisq(2 * (logLik(Poi2_select)[1]-logLik(NB2_select)[1]), 
       df = abs(NB2_select$df.residual - Poi2_select$df.residual), lower.tail=FALSE)


#Lack of fit testing
pchisq(q= sum(resid(Poi2_all,type='pearson')^2), df= summary(Poi1_all)$df.residual, lower.tail=FALSE)
pchisq(q= sum(resid(Poi2_int,type='pearson')^2), df= summary(Poi1_int)$df.residual, lower.tail=FALSE)
pchisq(q= sum(resid(Poi2_select,type='pearson')^2), df= summary(Poi1_select)$df.residual, lower.tail=FALSE)
pchisq(q= sum(resid(NB2_all,type='pearson')^2) , df= summary(NB2_all)$df.residual , lower.tail=FALSE)
pchisq(q= sum(resid(NB2_int,type='pearson')^2) , df= summary(NB2_int)$df.residual , lower.tail=FALSE)
pchisq(q= sum(resid(NB2_select,type='pearson')^2) , df= summary(NB2_select)$df.residual , lower.tail=FALSE)

#Information Criteria for All 3-5 models
AIC(Poi2_all);BIC(Poi2_all)
AIC(Poi2_int);BIC(Poi2_int)
AIC(Poi2_select);BIC(Poi2_select)
AIC(NB2_all);BIC(NB2_all)
AIC(NB2_int);BIC(NB2_int)
AIC(NB2_select);BIC(NB2_select)

#Intercept Only Poisson Model
Poi2_int <-glm(Freq ~ 1,data= GLM_Data,subset=Cat_HURDAT=='3-5', family=poisson)
#Intercept Only NB Model
NB2_int<-glm.nb(Freq ~ 1,data= GLM_Data,subset= Cat_HURDAT=='3-5', maxit = 1000)
Phi2_int<-summary(NB2_int)$theta

mu2<-mean(GLM_Data[GLM_Data$Cat_HURDAT=='3-5',]$Freq)

set.seed(1212)
test.data<-data.frame(Value=c(rpois(1000, mu2),rnegbin(1000, mu2, Phi2_int),GLM_Data[GLM_Data$Cat_HURDAT=='3-5',]$Freq,rep(NA, 1000-length(GLM_Data[GLM_Data$Cat_HURDAT=='3-5',]$Freq))), Type=c(rep('Pois. Sim.',1000),rep('NB Sim.',1000),rep('Empirical',1000)),check.names=F)
#pdf("Density35_new.pdf")
ggplot(data=test.data,aes(x= Value,y = ..density..,group=Type,fill=Type))+geom_bar(stat="bin",position= position_dodge(.5))+theme_bw()+ylab('Density')+scale_x_continuous(breaks = seq(0,11, by = 1),'# of Category 3-5 Tropical Cyclones')+ scale_fill_manual(values=c('grey0','grey50','grey75'))+theme(legend.position = c(0.7, 0.9), legend.justification = c(0, 1), panel.grid.major = element_blank(),axis.title=element_text(size=18),axis.text=element_text(size=16), legend.title=element_blank() ,legend.text=element_text(size=14), title =element_text(size=16))
#dev.off()


#ggplot(data=test.data,aes(x= Value,y = ..density..,group=Type,fill=Type))+geom_bar(stat="bin",position= position_dodge(.5))+ labs(title = "Empirical Density of Cat. 3-5 Tropical Cyclones\n vs. Poisson and NB Densities")+theme_bw()+ylab('Density')+scale_x_continuous(breaks = seq(0,11, by = 1),'# of Category 3-5 Tropical Cyclones')+ scale_fill_manual(values=c('#1b9e77','#d95f02','#7570b3'))+theme(legend.position = c(0.7, 0.9), legend.justification = c(0, 1), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.title=element_text(size=16),axis.text=element_text(size=14), legend.title=element_blank() ,legend.text=element_text(size=14), title =element_text(size=16))
#################################################################
#Setting Hyperparameters for Beta and phi
#################################################################

#Hyperparameters for Beta (mean and precision)
model.mat.x <-model.matrix(Poi1_all)
#model.mat.x.select <- model.matrix(Poi1_select)
model.mat.x1 <- model.matrix(Poi1_select)
model.mat.x2 <- model.matrix(Poi2_select)

#Largest Poi Model For Cat TS-2
hyper_Beta1_poi_all <-matrix(c(coef(Poi1_all),rep(1e-4,7)),nrow=2,byrow=T)

#Intercept only Poi Model For Cat TS-2
hyper_Beta1_poi_int <-matrix(c(coef(Poi1_int),rep(1e-4,1)),nrow=2,byrow=T)

#Selected Poi Model For Cat TS-2
hyper_Beta1_poi_selected <-matrix(c(coef(Poi1_select),rep(1e-4,5)),nrow=2,byrow=T)

#Largest NB Model For Cat TS-2
hyper_Beta1_nb_all <-matrix(c(coef(NB1_all),rep(1e-4,7)),nrow=2,byrow=T)
hyper_phi1_nb_all <-matrix(c(summary(NB1_all)$theta,rep(1,1)),nrow=2,byrow=T)

#Intercept only NB Model For Cat TS-2
hyper_Beta1_nb_int <-matrix(c(coef(NB1_int),rep(1e-4,1)),nrow=2,byrow=T)
hyper_phi1_nb_int <-matrix(c(summary(NB1_int)$theta,rep(1,1)),nrow=2,byrow=T)

#Selected Model For Cat TS-2
hyper_Beta1_nb_selected <-matrix(c(coef(NB1_select),rep(1e-4,5)),nrow=2,byrow=T)
hyper_phi1_nb_selected <-matrix(c(summary(NB1_select)$theta,rep(1,1)),nrow=2,byrow=T)

#Largest Poi Model For Cat 3-5
hyper_Beta2_poi_all <-matrix(c(coef(Poi2_all),rep(1e-4,7)),nrow=2,byrow=T)

#Intercept only Model For Cat 3-5
hyper_Beta2_poi_int <-matrix(c(coef(Poi2_int),rep(1e-4,1)),nrow=2,byrow=T)

#Selected Model For Cat 3-5
hyper_Beta2_poi_selected <-matrix(c(coef(Poi2_select),rep(1e-4,4)),nrow=2,byrow=T)


