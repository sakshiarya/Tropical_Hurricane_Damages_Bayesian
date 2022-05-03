#################### Bivariate logistic EVD ##################################################################
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
############################################################################################################
library(evd)

logpost_factory <- function(priorscale, data) function(par){
  n <- nrow(data)  
  ## mu1 nonstationary
  mu1 <- par[1] + par[2]*data$z1 + par[3]*data$z2 
  mu2 <- par[6] + par[7]*data$z1 + par[8]*data$z2
  
  ## t(x1)  
  term1_tx1 <- 1 + ((par[5]*(data$x1 - mu1))/(par[4]))
  tx1_exponent <- -(1/par[5])
  if(par[5]!= 0 & all(term1_tx1 > 0)){
    tx1 <- term1_tx1^tx1_exponent
  }
  else if(par[5] == 0 & all(term1_tx1 > 0)){
    tx1 <- exp(-(data$x1-mu1)/par[4])
  }
  else{
    tx1 <- 0
  }
  ## t(x2) 
  term1_tx2 <- 1 + ((par[10]*(data$x2 - mu2))/(par[9]))
  tx2_exponent <- -(1/par[10])
  if(par[10]!= 0 & all(term1_tx2 > 0)){
    tx2 <- term1_tx2^tx2_exponent
  }
  else if(par[10] == 0 & all(term1_tx2 > 0)){
    tx2 <- exp(-(data$x2-mu2)/par[9])
  }
  else{
    tx2 <- 0
  }
  if(par[4] > 0 & par[9] >0 & par[11] >0.05 & par[11] <=1 & all(tx1 != 0) & all(tx2 != 0)){
    loglik <- sum((par[5]+1)*log(tx1) + (par[10]+1)*log(tx2) + 
            ((1/par[11]) - 1)*log(tx1) + ((1/par[11])-1)*log(tx2) - 
              (tx1^(1/par[11]) + tx2^(1/par[11]))^par[11] + 
            (par[11]-2)*log(tx1^(1/par[11]) + tx2^(1/par[11])) + 
            log(((1-par[11])/par[11]) + (tx1^(1/par[11]) + tx2^(1/par[11]))^par[11])) +
            n*(-log(par[4]) - log(par[9]) - (par[1]^2)/(2*priorscale[1]) - 
           (par[2]^2/(2*priorscale[2]))- (par[3]^2/(2*priorscale[3])) -
           (priorscale[4] + 1)*log(par[4]) - (priorscale[5]/par[4]) -
           (par[6]^2)/(2*priorscale[8]) - (par[7]^2/(2*priorscale[9]))-
           (par[8]^2/(2*priorscale[10]))-
           (priorscale[11] + 1)*log(par[9]) - (priorscale[12]/par[9]) +
           log(ifelse(par[5] > priorscale[6] & par[5]< priorscale[7], 1/(priorscale[7]-priorscale[6]), 0)) + 
           log(ifelse(par[10] > priorscale[13] & par[10] < priorscale[14], 1/(priorscale[14] - priorscale[13]), 0)) + 
           log(ifelse(par[11]>= priorscale[15] & par[11]<= priorscale[16], 1/(priorscale[16] - priorscale[15]), 0)))
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


tri_var$x1 <- tri_var$logmaxWS
tri_var$x2 <- tri_var$logdamage
# normalize covariates
tri_var$z1 <- scale(tri_var$logminCP) 
tri_var$z2 <- scale(tri_var$avgLat)



## Scaling for hyperparameters
var_gam10 <- 100; var_gam11 <- 100; var_gam12 <- 100; alpha_sig1 <- 1; beta_sig1 <- 1;  a_xi1 <- -0.5; b_xi1 <- 0.5;
var_gam20 <- 1000; var_gam21 <- 1000; var_gam22 <- 1000; alpha_sig2 <- 1; beta_sig2 <- 3; a_xi2 <- -0.5; b_xi2 <- 0.5;
 a_r <- 0.05; b_r <- 1

priorscale <- c(var_gam10, var_gam11, var_gam12, alpha_sig1, beta_sig1,a_xi1, b_xi1,
                var_gam20, var_gam21, var_gam22,  alpha_sig2, beta_sig2, 
                a_xi2, b_xi2, a_r, b_r)

## step scales for MH-steps
step_scale <- c(0.03,rep(0.01,2),rep(0.005,2), 0.05, 0.04, 0.05,0.04, 0.02, 0.02)

## Fitting a parametric logistic bivariate model
mod1 <- fbvevd(cbind(tri_var$x1, tri_var$x2), nsloc1 = data.frame(tri_var$z1, tri_var$z2), nsloc2 = data.frame(tri_var$z1,tri_var$z2), model = "log")

pos_init <- mod1$estimate

logpost <- logpost_factory(data = tri_var, priorscale =  priorscale)

## Running the MH algorithm using the metrop function:
set.seed(608)
mo1 <- metrop(logpost, initial = pos_init, 1e6, scale = step_scale)
mo1$accept
results <- mo1$batch
post_mean <- apply(mo1$batch,2, mean)
post_sd <- apply(mo1$batch,2, sd)

mcse.results <- mcse.mat(mo1$batch)
save(results, file = "NS_MH_1e6_metrop_6August.RData")



## Plots diagnostics
load("NS_MH_1e6_metrop_6August.RData")
results_every100 <- results[seq(1, nrow(mo1$batch), 100),]
dim(results_every100)

colnames(results_every100) <- c("gamma10", "gamma11", "gamma12", "sigma1", "xi1", "gamma20",
                                "gamma21", "gamma22", "sigma2","xi2","r")
head(results_every100)

filefig <- "/Users/sakshi/Documents/School/Research/Lindsey_Dissertation2021/Chapter 1 - Wind Speed Pressure/BayesianBivEVDMCMC/NonStationaryBiv/FullNS/6August/"

 pdf(paste0(filefig,"TS_NSBivEVD_6August_1e6_1to5.pdf"))
 plot(ts(results_every100[,1:5]), main = "Time Series plots for MCMC output")
 dev.off()
 
 pdf(paste0(filefig,"TS_NSBivEVD_6August_1e6_6to11.pdf"))
 plot(ts(results_every100[,6:11]), main = "Time Series plots for MCMC output")
 dev.off() 
 
 pdf(paste0(filefig,"ACF_NSBivEVD_6August_1e6_1to3.pdf"))
 acf(results[,1:3], main = "ACF plot")
 dev.off()

 pdf(paste0(filefig,"ACF_NSBivEVD_6August_1e6_4to5.pdf"))
 acf(results[,4:5], main = "ACF plot")
 dev.off()
 
 pdf(paste0(filefig,"ACF_NSBivEVD_6August_1e6_9to11.pdf"))
 acf(results[,9:11], main = "ACF plot")
 dev.off()
 
 
xtable(cbind(post_mean, post_sd, mcse.results[,2], mod1$estimate, mod1$std.err), digits= 4)



