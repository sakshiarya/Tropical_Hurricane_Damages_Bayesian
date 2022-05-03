#################### Hierarchical NS EVD ##########################################################################
# This code is for fitting a hierarchical non-stationary EVD model with all three variables maxWS, damages, minCP as 
# univariate EVD models.

### logpost_factory:
#This function computes the log unnormalized posterior density (log-likelihood + log-prior)
### Input: 
# priorscale: Scaling for hyperparameters
# data: (x2, x1, z1, z2) i.e. (damages, maxWS, minCP, avgLat)

### Metrop function in mcmc package:
## Use to run the MH-algorithm
## Input:
# step_scale: step scales for MH-steps
# pos_init: the initial state of the Markov chain

## Author: Sakshi Arya (aryax010@umn.edu)
## Date : 4th August, 2021
####################################################################################################################



logpost_factory <- function(priorscale, data) function(par){
  n <- nrow(data)  
  ## mu1 nonstationary
  mu1 <- par[1] + par[2]*data$z2  
  mu2 <- par[5] + par[6]*data$z1 + par[7]*data$z2
  mu3 <- par[10] + par[11]*data$scalex1 + par[12]*data$z1 + par[13]*data$z2
  
  ## t(x1)  
  term1_tz1 <- 1 + ((par[4]*(data$z1 - mu1))/(par[3]))
  tz1_exponent <- -(1/par[4])
  if(par[4]!= 0 & all(term1_tz1 > 0)){
    tz1 <- term1_tz1^tz1_exponent
  }
  else if(par[4] == 0 & all(term1_tz1 > 0)){
    tz1 <- exp(-(data$z1-mu1)/par[3])
  }
  else{
    tz1 <- 0
  }
  ## t(x1) 
  term1_tx1 <- 1 + ((par[9]*(data$x1 - mu2))/(par[8]))
  tx1_exponent <- -(1/par[9])
  if(par[9]!= 0 & all(term1_tx1 > 0)){
    tx1 <- term1_tx1^tx1_exponent
  }
  else if(par[9] == 0 & all(term1_tx1 > 0)){
    tx1 <- exp(-(data$x1-mu2)/par[8])
  }
  else{
    tx1 <- 0
  }
  ## t(x2) 
  term1_tx2 <- 1 + ((par[15]*(data$x2 - mu3))/(par[14]))
  tx2_exponent <- -(1/par[15])
  if(par[15]!= 0 & all(term1_tx2 > 0)){
    tx2 <- term1_tx2^tx2_exponent
  }
  else if(par[15] == 0 & all(term1_tx2 > 0)){
    tx2 <- exp(-(data$x2-mu3)/par[14])
  }
  else{
    tx2 <- 0
  }
  loglik_EVD <- function(y, sigma, xi){
     -log(sigma) + (xi + 1)*log(y) - y
  }

  if(par[3] > 0 & par[8] >0 & par[14] > 0 & all(tx1 != 0) & all(tx2 != 0)){
    loglik <- sum(loglik_EVD(y = tz1, sigma = par[3], xi = par[4]) + 
                    loglik_EVD(y = tx1, sigma = par[8], xi = par[9]) + 
                    loglik_EVD(y = tx2, sigma = par[14], xi = par[15])) +
      n*( - (par[1]^2/(2*priorscale[1]))- (par[2]^2/(2*priorscale[2])) -
            (par[5]^2/(2*priorscale[7])) - (par[6]^2/(2*priorscale[8])) -
            (par[7]^2/(2*priorscale[9])) - (par[10]^2/(2*priorscale[14])) - 
            (par[11]^2/(2*priorscale[15])) - (par[12]^2/(2*priorscale[16])) -
           (par[13]^2/(2*priorscale[17])) - 
           (priorscale[3] + 1)*log(par[3]) - (priorscale[4]/par[3]) -
           (priorscale[10] + 1)*log(par[8]) - (priorscale[11]/par[8]) -
           (priorscale[18] + 1)*log(par[14]) - (priorscale[19]/par[14]) +
           log(ifelse(par[4] > priorscale[5] & par[4]< priorscale[6], 1/(priorscale[6]-priorscale[5]), 0)) + 
           log(ifelse(par[9] > priorscale[12] & par[9] < priorscale[13], 1/(priorscale[13] - priorscale[12]), 0)) + 
           log(ifelse(par[15]>= priorscale[20] & par[15]<= priorscale[21], 1/(priorscale[21] - priorscale[20]), 0)))
  }else{
    loglik = -Inf
  }
  return(loglik)
}


## Reading the data:
storms<-read.csv('Categorized_Storm_1960_2019.csv')

names(storms)
storm_damage<-storms[storms$minCP< Inf & storms$CURRENT.DAMAGE.2021>0 & !is.na(storms$CURRENT.DAMAGE.2021),]

tri_var<-storm_damage[,c('minCP','maxWS','CURRENT.DAMAGE.2021')]
tri_var$logdamage<-log(tri_var$CURRENT.DAMAGE.2021)
tri_var$logmaxWS<-log(tri_var$maxWS)

tri_var$logminCP<-log(1013-tri_var$minCP)
tri_var$avgLat <- storm_damage$avgLat

densityPlot(tri_var$logdamage)

tri_var$x1 <- tri_var$logmaxWS
tri_var$scalex1 <- scale(tri_var$logmaxWS)
tri_var$x2 <- tri_var$logdamage
tri_var$z1 <- scale(tri_var$logminCP)
tri_var$z2 <- scale(tri_var$avgLat)

# Frequentist models
m_z1 <- fevd(z1, data = tri_var,location.fun=~z2, scale.fun=~1, shape.fun=~1,type='GEV')
m_z1$results$par
m_x1 <- fevd(x1, data = tri_var,location.fun=~ z1 + z2, scale.fun=~1, shape.fun=~1,type='GEV')
m_x1$results$par
m_x2 <- fevd(x2, data = tri_var,location.fun=~ scalex1 + z1 + z2, scale.fun=~1, shape.fun=~1,type='GEV')
m_x2$results$par

freq_se <- c(summary(m_z1)$se.theta, summary(m_x1)$se.theta, summary(m_x2)$se.theta)


pos_init <- c(m_z1$results$par, m_x1$results$par, m_x2$results$par)
names(pos_init)
xtable(as.data.frame(pos_init), digits = 4)



## Scaling for hyperparameters
var_alpha0 <- 100; var_alpha1 <- 100; alpha_sigz1 <- 1; beta_sigz1 <- 1; a_z1 <- -1; b_z1 <- 1
var_beta0 <- 100; var_beta1 <- 100;  var_beta2 <- 100; alpha_sigx1 <- 1; beta_sigx1 <- 1;  a_x1 <- -0.5; b_x1 <- 0.5
var_gam0 <- 10^3; var_gam1 <- 10^3; var_gam2 <- 10^3; var_gam3 <- 10^3; alpha_sigx2 <- 2; beta_sigx2 <- 3; a_x2 <- -0.5; b_x2 <- 0.5


priorscale <- c(var_alpha0, var_alpha1, alpha_sigz1, beta_sigz1, a_z1, b_z1, 
                var_beta0, var_beta1, var_beta2, alpha_sigx1, beta_sigx1, a_x1, b_x1,
                var_gam0, var_gam1, var_gam2, var_gam3, alpha_sigx2, beta_sigx2, a_x2, b_x2)

step_scale <- c(rep(0.01, 9), rep(0.09,4), rep(0.02,2))


logpost <- logpost_factory(data = tri_var, priorscale =  priorscale)
set.seed(1708)
mo1 <- metrop(logpost, initial = pos_init, 1e6, scale = step_scale)
mo1$accept
results <- mo1$batch


save(results, file = "hierNS_1e6_6August.RData")
colnames(results) <- c("alpha0", "alpha1","sigma_z1","xi_z1", "beta0","beta1","beta2",
                       "sigma_x1", "xi_x1", "gamma0", "gamma1", "gamma2","gamma3",
                       "sigma_x2", "xi_x2")

results <- mo1$batch


load("hierNS_1e6_6August.RData")
post_means <- apply(results, 2, mean)
post_sd <- apply(results, 2, sd)
xtable(cbind(post_means, post_sd), digits = 4)
results_every100 <- results[seq(1, nrow(mo1$batch), 100),]
dim(results_every100)

filefig <- "/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/BayesianBivEVDMCMC/NonStationaryBiv/HierarchicalNS/6Aug/"
pdf(paste0(filefig,"TS_Hier_6August_1to4.pdf"))
plot.ts(results_every100[,1:4], main = "Time series plots for the MCMC output")
dev.off()
pdf(paste0(filefig,"TS_Hier_6August_5to9.pdf"))
plot.ts(results_every100[,5:9], main = "Time series plots for the MCMC output")
dev.off()
pdf(paste0(filefig,"TS_Hier_6August_10to15.pdf"))
plot.ts(results_every100[,10:15], main = "Time series plots for the MCMC output")
dev.off()

pdf(paste0(filefig,"ACF_Hier_6August_1to4.pdf"))
acf(results[,1:4], main = "ACF plot for the MCMC output")
dev.off()
pdf(paste0(filefig,"ACF_Hier_6August_5to9.pdf"))
acf(results[,5:9], main = "ACF plot for the MCMC output")
dev.off()
pdf(paste0(filefig,"ACF_Hier_6August_10to13.pdf"))
acf(results[,10:13], main = "ACF plot for the MCMC output")
dev.off()
pdf(paste0(filefig,"ACF_Hier_6August_14to15.pdf"))
acf(results[,14:15], main = "ACF plot for the MCMC output")
dev.off()

apply(mo1$batch, 2, mean)
pos_init
plot.ts(mo1$batch[,1:9])
plot.ts(mo1$batch[,10:15])
acf(mo1$batch[,10:14])

