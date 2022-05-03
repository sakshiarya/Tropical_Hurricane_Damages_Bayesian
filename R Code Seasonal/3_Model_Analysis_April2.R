#################################################################
#################################################################
#Title: Bayesian Hierarchical Models diagnostics
#Author: Lindsey Dietz (diet0146@umn.edu)
#Objective: Create models and run MCMC
#Created: 2/25/16
# Last Updated: 4/2/2021 (by Sakshi Arya)
#################################################################
#################################################################

#Run if needed
filep<-'/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 2 - Bayesian Damage/'
#source(paste(filep,'0_Data_Processing.r',sep=''))
#source(paste(filep,'1a_Hyperpar_Beta.r',sep=''))
#source(paste(filep,'1b_Hyperpar_Theta.r',sep=''))
#source(paste(filep,'1c_Hyperpar_Mu_Sigma.r',sep=''))
#source(paste(filep,'2_Models.r',sep=''))

#Loads all necessary libraries
library(mcmcse)
library(coda)
library(dclone)
library(mcmcplots)
library(mvtnorm)
library(xtable)

#Load results of MCMC
load(paste(filep,'MCMC Chains/March 31/Full Bayesian/Select_covs/clones.MCMC.final.FB.select.April5_withSST.rda',sep=''))
load(paste(filep,'MCMC Chains/March 31/Full Bayesian/Select_covs/jags.MCMC.final.FB.select.April5_withSST.rda',sep=''))

#################################################################
#Dclone Model
#################################################################
summary(clones.MCMC.final.FB.select.April5_withSST)
dcdiag(clones.MCMC.final.EB.all.April1)
confint(clones.MCMC.final.EB.all.April1)
coef(clones.MCMC.final.EB.all.April1)
dcsd(clones.MCMC.final.EB.all.April1)

xtable(summary(clones.MCMC.final.FB.select.April5_withSST)[1]$statistics, digits = c(3,3,3,3,-3,-3,3))

par(mfrow=c(3,7))
densplot(clones.MCMC.final.FB.select.April5_withSST, xlab='')
autocorr.plot(clones.MCMC.final.FB.select.April5_withSST)
traceplot(clones.MCMC.final.FB.select.April5_withSST)

#################################################################
#Jags Model
#################################################################
jags.MCMC.final.EB.all.April1$BUGSoutput$summary
xtable(jags.MCMC.final.EB.all.April1$BUGSoutput$summary)

jags.posterior <- as.mcmc(x= jags.MCMC.final.EB.all.April1)
summary(jags.posterior)

multi.se.EB.all <- mcse.mat(rbind(jags.posterior[[1]][,c(-6)],jags.posterior[[2]][,c(-6)],jags.posterior[[3]][,c(-6)]), method = "bm", g = NULL)
xtable(multi.se.all, digits = c(-3,-3,-3))

xtable(summary(jags.posterior)[1]$statistics, digits = c(3,3,3,-3,-3))
cat <- c(rep(1,6),rep(2,6))
ind_est <- c(seq(1,6,by=1),seq(1,6,by=1))
for(i in 1:12){
  par(mfrow=c(1,1))
  #pdf(paste0("beta",cat[i],ind_est[i],"_acf_EB_full.pdf"),height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main = substitute(paste('Autocorrelation Function of ', beta[m][n]), list(m=cat[i], n =ind_est[i])), cex.lab = 1.5, cex.main=1.5)
  #dev.off()
}

for(i in 13:14){
  par(mfrow = c(1,1))
  #pdf(paste0("theta",i-12,"_acf_EB_full.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',theta[m]), list(m=c(1,2)[i-12])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 16:17){
  par(mfrow = c(1,1))
  #pdf(paste0("mu",i-15,"_acf_EB_full.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',mu[m]), list(m=c(1,2)[i-15])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 18:19){
  par(mfrow = c(1,1))
  #pdf(paste0("sigma",i-17,"_acf_EB_full.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',sigma[m]), list(m=c(1,2)[i-17])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

par(mfrow=c(1,1))
#pdf("r_acf_EB_full.pdf", height = 6, width = 8)
autocorr.plot(x=jags.posterior[[1]][,c(20)], auto.layout = F, main=expression(paste('Autocorrelation Function of ',r)), cex.lab=1.5, cex.main=1.5)
#dev.off()


#################################################################
#################################################################
#Trace and Density Plots for Each Parameter in TS-2
#################################################################
#################################################################
#Beta01
for(i in 1:12){
  #pdf(paste0("beta",cat[i],ind_est[i],"_EB_full.pdf"),height = 5, width = 10)
  par(mfrow=c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main = substitute(paste('Trace of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main = substitute(paste('Density of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 13:14){
  #pdf(paste0("theta",i-12,"_EB_full.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',theta[m]), list(m=c(1,2)[i-12])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',theta[m]), list(m=c(1,2)[i-12])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 16:17){
  #pdf(paste0("mu",i-15,"_EB_full.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',mu[m]), list(m=c(1,2)[i-15])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',mu[m]), list(m=c(1,2)[i-15])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 18:19){
  #pdf(paste0("sigma",i-17,"_EB_full.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',sigma[m]), list(m=c(1,2)[i-17])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',sigma[m]), list(m=c(1,2)[i-17])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

#r (overdispersion parameter)
#pdf("r_EB_full.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
traceplot(jags.posterior[[1]][,c(20)],main=expression(paste('Trace of ',phi)), xlab='Iteration', ylab='Estimate',cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
densplot(jags.posterior[[1]][,c(20)],main=expression(paste('Density of ',phi)), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
#dev.off()


multi.cov.all <- mcse.multi(rbind(jags.posterior[[1]][,c(-6)],jags.posterior[[2]][,c(-6)],jags.posterior[[3]][,c(-6)]), blather = TRUE)
multi.se.all <- mcse.mat(rbind(jags.posterior[[1]][,c(-6)],jags.posterior[[2]][,c(-6)],jags.posterior[[3]][,c(-6)]), method = "bm", g = NULL)
diag(multi.cov.all$cov/multi.cov.all$size)
multiESS(rbind(jags.posterior[[1]][,c(-6)],jags.posterior[[2]][,c(-6)],jags.posterior[[3]][,c(-6)]), covmat= multi.cov.all)


##############################################################################
###############################################################################
### Empirical Bayes With selected model covariates in preliminary analysis ###
###############################################################################

#Load results of MCMC
load(paste(filep,'MCMC Chains/March 31/Empirical Bayesian/Selected_covs/clones.MCMC.final.EB.select.April5.SST.rda',sep=''))
load(paste(filep,'MCMC Chains/March 31/Empirical Bayesian/Selected_covs/jags.MCMC.final.EB.select.April5.SST.rda',sep=''))

#################################################################
#Dclone Model
#################################################################
summary(clones.MCMC.final.EB.select.April5.SST)
dcdiag(clones.MCMC.final.EB.select.April5.SST)
confint(clones.MCMC.final.EB.select.April5.SST)
coef(clones.MCMC.final.EB.select.April5.SST)
dcsd(clones.MCMC.final.EB.select.April5.SST)

xtable(summary(clones.MCMC.final.EB.select.April5.SST)[1]$statistics, digits = c(3,3,3,3,-3,-3,3))

par(mfrow=c(3,7))
densplot(clones.MCMC.final.EB.select.April5.SST, xlab='')
autocorr.plot(clones.MCMC.final.EB.select.April5.SST)
traceplot(clones.MCMC.final.EB.select.April5.SST)

#################################################################
#Jags Model
#################################################################
jags.MCMC.final.EB.select.April5.SST$BUGSoutput$summary
xtable(jags.MCMC.final.EB.select.April5.SST$BUGSoutput$summary)

jags.posterior <- as.mcmc(x= jags.MCMC.final.EB.select.April5.SST)
summary(jags.posterior)
multi.se.EB.select <- mcse.mat(rbind(jags.posterior[[1]][,c(-6)],jags.posterior[[2]][,c(-6)],jags.posterior[[3]][,c(-6)]), method = "bm", g = NULL)

xtable(summary(jags.posterior)[1]$statistics, digits = c(3,3,3,-3,-3))
cat <- c(rep(1,4),rep(2,3))
ind_est <- c(seq(1,4,by=1),seq(1,3,by=1))
for(i in 1:7){
  #pdf(paste0("beta",cat[i],ind_est[i],"_acf_EB_select.pdf"),height = 6, width = 8)
  par(mfrow=c(1,1))
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main = substitute(paste('Autocorrelation Function of ', beta[m][n]), list(m=cat[i], n =ind_est[i])), cex.lab = 1.5, cex.main=1.5)
  #dev.off()
}

for(i in 8:9){
  par(mfrow = c(1,1))
  #pdf(paste0("theta",i-7,"_acf_EB_select.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',theta[m]), list(m=c(1,2)[i-7])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 11:12){
  par(mfrow = c(1,1))
  #pdf(paste0("mu",i-10,"_acf_EB_select.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',mu[m]), list(m=c(1,2)[i-10])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 13:14){
  par(mfrow = c(1,1))
  #pdf(paste0("sigma",i-12,"_acf_EB_select.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',sigma[m]), list(m=c(1,2)[i-12])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

par(mfrow=c(1,1))
#pdf("r_acf_EB_select.pdf", height = 6, width = 8)
autocorr.plot(x=jags.posterior[[1]][,c(15)], auto.layout = F, main=expression(paste('Autocorrelation Function of ',r)), cex.lab=1.5, cex.main=1.5)
#dev.off()


#################################################################
#################################################################
#Trace and Density Plots for Each Parameter in TS-2 and Cat 3-5
#################################################################
#################################################################
#Beta01
for(i in 1:7){
  #pdf(paste0("beta",cat[i],ind_est[i],"_EB_select.pdf"),height = 5, width = 10)
  par(mfrow=c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main = substitute(paste('Trace of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main = substitute(paste('Density of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 8:9){
  #pdf(paste0("theta",i-7,"_EB_select.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',theta[m]), list(m=c(1,2)[i-7])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',theta[m]), list(m=c(1,2)[i-7])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 11:12){
  #pdf(paste0("mu",i-10,"_EB_select.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',mu[m]), list(m=c(1,2)[i-10])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',mu[m]), list(m=c(1,2)[i-10])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 13:14){
  #pdf(paste0("sigma",i-12,"_EB_select.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',sigma[m]), list(m=c(1,2)[i-12])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',sigma[m]), list(m=c(1,2)[i-12])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

#r
#pdf("r_EB_select.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
traceplot(jags.posterior[[1]][,c(15)],main=expression(paste('Trace of ',r)), xlab='Iteration', ylab='Estimate',cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
densplot(jags.posterior[[1]][,c(15)],main=expression(paste('Density of ',r)), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
#dev.off()

##############################################################################
###############################################################################
###################### Fully Bayes with all covariates  ######################
###############################################################################
#rm(list = ls())
#Load results of MCMC
load(paste(filep,'MCMC Chains/Full Bayes/March 31/Full Bayesian/All_covs/clones.MCMC.final.FB.all.April2.rda',sep=''))
load(paste(filep,'MCMC Chains/Full Bayes/March 31/Full Bayesian/All_covs/jags.MCMC.final.FB.all.April2.rda',sep=''))

#################################################################
#Dclone Model
#################################################################
summary(clones.MCMC.final.FB.all.April2) # The EB in the name is mis-spelled it should be FB
dcdiag(clones.MCMC.final.FB.all.April2)
confint(clones.MCMC.final.FB.all.April2)
coef(clones.MCMC.final.FB.all.April2)
dcsd(clones.MCMC.final.FB.all.April2)

xtable(summary(clones.MCMC.final.FB.all.April2)[1]$statistics, digits = c(3,3,3,3,-3,-3,3))
xtable(summary(clones.MCMC.final.FB.all.April2)[1]$statistics[,1:2], digits = c(3,3,3))

par(mfrow=c(3,7))
densplot(clones.MCMC.final.FB.all.April2, xlab='')
autocorr.plot(clones.MCMC.final.FB.all.April2)
traceplot(clones.MCMC.final.FB.all.April2)

#################################################################
#Jags Model
#################################################################
jags.MCMC.final.FB.all.April2$BUGSoutput$summary
xtable(jags.MCMC.final.FB.all.April2$BUGSoutput$summary)

jags.posterior <- as.mcmc(x= jags.MCMC.final.FB.all.April2)
summary(jags.posterior)
multi.se.FB.all <- mcse.mat(rbind(jags.posterior[[1]][,c(-6)],jags.posterior[[2]][,c(-6)],jags.posterior[[3]][,c(-6)]), method = "bm", g = NULL)

xtable(summary(jags.posterior)[1]$statistics[,1:2], digits = c(3,3,3))
cat <- c(rep(1,6),rep(2,6))
ind_est <- c(seq(1,6,by=1),seq(1,6,by=1))
seq_100 <- seq(1,99900,by =100)
for(i in 1:12){
  par(mfrow=c(1,1))
  #pdf(paste0("beta",cat[i],ind_est[i],"_acf_BayesFfull.pdf"),height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][seq_100,c(i)], auto.layout = F, main = substitute(paste('Autocorrelation Function of ', beta[m][n]), list(m=cat[i], n =ind_est[i]-1)), cex.lab = 1.5, cex.main=1.5)
  #dev.off()
}

for(i in 13:14){
  par(mfrow = c(1,1))
  #pdf(paste0("theta",i-12,"_acf_BayesFfull.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][seq_100,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',theta[m]), list(m=c(1,2)[i-12])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 16:17){
  par(mfrow = c(1,1))
  #pdf(paste0("mu",i-15,"_acf_BayesFfull.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',mu[m]), list(m=c(1,2)[i-15])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 18:19){
  par(mfrow = c(1,1))
  #pdf(paste0("sigma",i-17,"_acf_BayesFfull.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',sigma[m]), list(m=c(1,2)[i-17])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

par(mfrow=c(1,1))
#pdf("r_acf_BayesFfull.pdf", height = 6, width = 8)
autocorr.plot(x=jags.posterior[[1]][,c(20)], auto.layout = F, main=expression(paste('Autocorrelation Function of ',r)), cex.lab=1.5, cex.main=1.5)
#dev.off()


#################################################################
#################################################################
#Trace and Density Plots for Each Parameter in TS-2 and 3-5
#################################################################
#################################################################
#Beta01
for(i in 1:12){
  #pdf(paste0("beta",cat[i],ind_est[i],"_BayesFfull.pdf"),height = 5, width = 10)
  par(mfrow=c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main = substitute(paste('Trace of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main = substitute(paste('Density of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 13:14){
  #pdf(paste0("theta",i-12,"_BayesFfull.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',theta[m]), list(m=c(1,2)[i-12])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',theta[m]), list(m=c(1,2)[i-12])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 16:17){
  #pdf(paste0("mu",i-15,"_BayesFfull.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',mu[m]), list(m=c(1,2)[i-15])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',mu[m]), list(m=c(1,2)[i-15])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 18:19){
  #pdf(paste0("sigma",i-17,"_BayesFfull.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',sigma[m]), list(m=c(1,2)[i-17])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',sigma[m]), list(m=c(1,2)[i-17])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

#r
#pdf("r_BayesFfull.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
traceplot(jags.posterior[[1]][,c(20)],main=expression(paste('Trace of ',r)), xlab='Iteration', ylab='Estimate',cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
densplot(jags.posterior[[1]][,c(20)],main=expression(paste('Density of ',r)), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
#dev.off()



##############################################################################
###############################################################################
### Fully Bayes With selected model covariates in preliminary analysis ###
###############################################################################

#Load results of MCMC
load(paste(filep,'MCMC Chains/March 31/Full Bayesian/Select_covs/clones.MCMC.final.FB.select.April5_withSST.rda',sep=''))
load(paste(filep,'MCMC Chains/March 31/Full Bayesian/Select_covs/jags.MCMC.final.FB.select.April5_withSST.rda',sep=''))

xtable(summary(clones.MCMC.final.FB.select.April5_withSST)[1]$statistics, digits = c(3,3,3,3,-3,-3,3))

#################################################################
#Jags Model
#################################################################
jags.MCMC.final.FB.select.April5_withSST$BUGSoutput$summary
xtable(jags.MCMC.final.FB.select.April5_withSST$BUGSoutput$summary)

jags.posterior <- as.mcmc(x= jags.MCMC.final.FB.select.April5_withSST)
summary(jags.posterior)
multi.se.FB.select <- mcse.mat(rbind(jags.posterior[[1]][,c(-6)],jags.posterior[[2]][,c(-6)],jags.posterior[[3]][,c(-6)]), method = "bm", g = NULL)

xtable(summary(jags.posterior)[1]$statistics, digits = c(3,3,3,-3,-3))
cat <- c(rep(1,4),rep(2,3))
ind_est <- c(seq(1,4,by=1),seq(1,3,by=1))
for(i in 1:7){
  #pdf(paste0("beta",cat[i],ind_est[i],"_acf_BayesF_select.pdf"),height = 6, width = 8)
  par(mfrow=c(1,1))
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main = substitute(paste('Autocorrelation Function of ', beta[m][n]), list(m=cat[i], n =ind_est[i])), cex.lab = 1.5, cex.main=1.5)
  #dev.off()
}

for(i in 8:9){
  par(mfrow = c(1,1))
  #pdf(paste0("theta",i-7,"_acf_BayesF_select.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',theta[m]), list(m=c(1,2)[i-7])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 11:12){
  par(mfrow = c(1,1))
  #pdf(paste0("mu",i-10,"_acf_BayesF_select.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',mu[m]), list(m=c(1,2)[i-10])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

for(i in 13:14){
  par(mfrow = c(1,1))
  #pdf(paste0("sigma",i-12,"_acf_BayesF_select.pdf"), height = 6, width = 8)
  autocorr.plot(x=jags.posterior[[1]][,c(i)], auto.layout = F, main=substitute(paste('Autocorrelation Function of ',sigma[m]), list(m=c(1,2)[i-12])), cex.lab=1.5, cex.main=1.5)
  #dev.off()
}

par(mfrow=c(1,1))
#pdf("r_acf_BayesF_select.pdf", height = 6, width = 8)
autocorr.plot(x=jags.posterior[[1]][,c(15)], auto.layout = F, main=expression(paste('Autocorrelation Function of ',r)), cex.lab=1.5, cex.main=1.5)
#dev.off()


#################################################################
#################################################################
#Trace and Density Plots for Each Parameter in TS-2 and Cat 3-5
#################################################################
#################################################################
#Beta01
for(i in 1:7){
  #pdf(paste0("beta",cat[i],ind_est[i],"_BayesF_select.pdf"),height = 5, width = 10)
  par(mfrow=c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main = substitute(paste('Trace of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main = substitute(paste('Density of ', beta[m][n]), list(m=cat[i], n =ind_est[i])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 8:9){
  #pdf(paste0("theta",i-7,"_BayesF_select.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',theta[m]), list(m=c(1,2)[i-7])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',theta[m]), list(m=c(1,2)[i-7])),  xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 11:12){
  #pdf(paste0("mu",i-10,"_BayesF_select.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',mu[m]), list(m=c(1,2)[i-10])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',mu[m]), list(m=c(1,2)[i-10])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

for(i in 13:14){
  #pdf(paste0("sigma",i-12,"_BayesF_select.pdf"), height = 5, width = 10)
  par(mfrow = c(1,2))
  traceplot(x=jags.posterior[[1]][,c(i)], main=substitute(paste('Trace of ',sigma[m]), list(m=c(1,2)[i-12])), xlab='Iteration', ylab='Estimate', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  densplot(jags.posterior[[1]][,c(i)], main=substitute(paste('Density of ',sigma[m]), list(m=c(1,2)[i-12])), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
  #dev.off()
}

#Phi
#pdf("r_BayesF_select.pdf", height = 5, width = 10)
par(mfrow = c(1,2))
traceplot(jags.posterior[[1]][,c(15)],main=expression(paste('Trace of ',r)), xlab='Iteration', ylab='Estimate',cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
densplot(jags.posterior[[1]][,c(15)],main=expression(paste('Density of ',r)), xlab='Estimate', ylab='Density', cex = 0.5, cex.axis = 1, cex.lab = 1.25, cex.main=1.25)
#dev.off()



# Multivariate SE:
multi_SE <- cbind(multi.se.EB.all[,2], multi.se.FB.all[,2])

cbind(multi.se.EB.select[,2], multi.se.FB.select[,2])
