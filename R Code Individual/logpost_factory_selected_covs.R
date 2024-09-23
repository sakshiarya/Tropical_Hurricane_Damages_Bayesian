### Function for Posterior in a hierarchical GEV/trivariate/GEV with log-Normal damages model with  damages, windspeed and pressure as variables in the hierarchy.
## Input: # priorscale: hyperparameter choices,
# data: to be inputted
# covs: covariates common to each layer of hierarchy
# model: either "Hier_GEV" or "Tri_GEV" or "HalfNormalDamagesEVD" (depending on the model being used)
## Output: # unnormalized posterior from the log-likelihood and priors

logpost_factory <- function(priorscale, data, covs) function(par, model){
  n <- nrow(data)  
  if(model == "Hier_GEV"){
  ## nonstationary location parameters
  mu1 <- par[1] + rowSums(t(t(data[,covs[[1]]])*par[2:(length(covs[[1]])+1)])) 
  mu2 <- par[length(covs[[1]]) + 2] + rowSums(t(t(data[,covs[[2]]])*par[(length(covs[[1]]) +length(covs[[2]]) + 2)]))
  mu3 <- par[(length(covs[[1]]) +length(covs[[2]]) + 3)] +  rowSums(t(t(data[,covs[[3]]])*par[(length(covs[[1]]) +length(covs[[2]]) + 4):(length(covs[[1]]) +length(covs[[2]]) + 3 + length(covs[[3]]))]))
  
  d <- (length(covs[[1]]) +length(covs[[2]]) + 3 + length(covs[[3]]))
  ## t(z1)  
  term1_tz1 <- 1 + ((par[d+1]*(data$z1 - mu1))/(par[d+2]))
  tz1_exponent <- -(1/par[d+1])
  if(par[d+1]!= 0 & all(term1_tz1 > 0)){
    tz1 <- term1_tz1^tz1_exponent
  }
  else if(par[d+1] == 0 & all(term1_tz1 > 0)){
    tz1 <- exp(-(data$z1-mu1)/par[d+2])
  }
  else{
    tz1 <- 0
  }
  ## t(x1) 
  term1_tx1 <- 1 + ((par[d+3]*(data$x1 - mu2))/(par[d+4]))
  tx1_exponent <- -(1/par[d+3])
  if(par[d+3]!= 0 & all(term1_tx1 > 0)){
    tx1 <- term1_tx1^tx1_exponent
  }
  else if(par[d+3] == 0 & all(term1_tx1 > 0)){
    tx1 <- exp(-(data$x1-mu2)/par[d+4])
  }
  else{
    tx1 <- 0
  }
  ## t(x2) 
  term1_tx2 <- 1 + ((par[d+5]*(data$x2 - mu3))/(par[d+6]))
  tx2_exponent <- -(1/par[d+5])
  if(par[d+5]!= 0 & all(term1_tx2 > 0)){
    tx2 <- term1_tx2^tx2_exponent
  }
  else if(par[d+5] == 0 & all(term1_tx2 > 0)){
    tx2 <- exp(-(data$x2-mu3)/par[d+6])
  }
  else{
    tx2 <- 0
  }
  loglik_EVD <- function(y, sigma, xi){
    -log(sigma) + (xi + 1)*log(y) - y
  }
  
  if(par[d+2] > 0 & par[d+4] >0 & par[d+6] > 0 & all(tz1 != 0) & all(tx1 != 0) & all(tx2 != 0)){
    loglik <- sum(loglik_EVD(y = tz1, sigma = par[d+2], xi = par[d+1]) + 
                    loglik_EVD(y = tx1, sigma = par[d+4], xi = par[d+3]) + 
                    loglik_EVD(y = tx2, sigma = par[d+6], xi = par[d+5])) +
      n*( - sum((par[1:(d)]^2)/(2*(priorscale[1:(d)])))- 
            (priorscale[d + 1] + 1)*log(par[d+2]) - (priorscale[d + 2]/par[d+2]) -
            (priorscale[(d) + 3] + 1)*log(par[d+4]) - (priorscale[(d) + 4]/par[d+4]) -
            (priorscale[(d) + 5] + 1)*log(par[d+6]) - (priorscale[(d) + 6]/par[d+6]) +
            log(ifelse(par[d+1] > priorscale[(d) + 7] & par[d+1]< priorscale[(d) + 8], 1/(priorscale[(d) + 8]-priorscale[(d) + 7]), 0)) + 
            log(ifelse(par[d+3] > priorscale[(d) + 9] & par[d+3] < priorscale[(d) + 10], 1/(priorscale[(d) + 10] - priorscale[(d) + 9]), 0)) + 
            log(ifelse(par[d+5]> priorscale[(d) + 11] & par[d+5]< priorscale[(d) + 12], 1/(priorscale[(d) + 12] - priorscale[(d) + 11]), 0)))
    
  }else{
    loglik = -Inf
  }}
  else if(model == "Tri_GEV"){
    ## nonstationary location parameters
    mu1 <- par[1] + rowSums(t(t(data[,covs[[1]]])*par[2:(length(covs[[1]])+1)])) 
    mu2 <- par[length(covs[[1]]) + 2] + rowSums(t(t(data[,covs[[2]]])*par[(length(covs[[1]]) +length(covs[[2]]) + 2)]))
    mu3 <- par[(length(covs[[1]]) +length(covs[[2]]) + 3)] +  rowSums(t(t(data[,covs[[3]]])*par[(length(covs[[1]]) +length(covs[[2]]) + 4):(length(covs[[1]]) +length(covs[[2]]) + 3 + length(covs[[3]]))]))
    
    d <- (length(covs[[1]]) +length(covs[[2]]) + 3 + length(covs[[3]]))
    ## t(z1)  
    term1_tz1 <- 1 + ((par[d+1]*(data$z1 - mu1))/(par[d+2]))
    tz1_exponent <- -(1/par[d+1])
    if(par[d+1]!= 0 & all(term1_tz1 > 0)){
      tz1 <- term1_tz1^tz1_exponent
    }
    else if(par[d+1] == 0 & all(term1_tz1 > 0)){
      tz1 <- exp(-(data$z1-mu1)/par[d+2])
    }
    else{
      tz1 <- 0
    }
    ## t(x1) 
    term1_tx1 <- 1 + ((par[d+3]*(data$x1 - mu2))/(par[d+4]))
    tx1_exponent <- -(1/par[d+3])
    if(par[d+3]!= 0 & all(term1_tx1 > 0)){
      tx1 <- term1_tx1^tx1_exponent
    }
    else if(par[d+3] == 0 & all(term1_tx1 > 0)){
      tx1 <- exp(-(data$x1-mu2)/par[d+4])
    }
    else{
      tx1 <- 0
    }
    ## t(x2) 
    term1_tx2 <- 1 + ((par[d+5]*(data$x2 - mu3))/(par[d+6]))
    tx2_exponent <- -(1/par[d+5])
    if(par[d+5]!= 0 & all(term1_tx2 > 0)){
      tx2 <- term1_tx2^tx2_exponent
    }
    else if(par[d+5] == 0 & all(term1_tx2 > 0)){
      tx2 <- exp(-(data$x2-mu3)/par[d+6])
    }
    else{
      tx2 <- 0
    }
    if(par[d+2] > 0 & par[d+4] >0 & par[d+6] > 0 & all(tz1 != 0) & all(tx1 != 0) & all(tx2 != 0) & par[d+7] > 0.05 & par[d+7] < 1){
      loglik <- sum((par[d+1]+1)*log(tz1) + (par[d+3]+1)*log(tx1) + (par[d+5]+1)*log(tx2) +
                      ((1/par[d+7]) - 1)*log(tz1) + ((1/par[d+7])-1)*log(tx1) + ((1/par[d+7])-1)*log(tx2) - 
                      (tx1^(1/par[d+7]) + tx2^(1/par[d+7]) + tz1^(1/par[d+7]))^par[d+7] + 
                      (par[(d)+7]-3)*log(tz1^(1/par[(d)+7]) + tx1^(1/par[(d)+7]) + tx2^(1/par[(d)+7])) + 
                      log((((1-par[(d)+7])*(2-par[(d)+7]))/par[(d)+7]) - ((1-par[(d)+7])/par[(d)+7])*(tz1^(1/par[(d)+7]) +tx1^(1/par[(d)+7]) + tx2^(1/par[(d)+7]))^par[(d)+7] + 
                            (tz1^(1/par[(d)+7]) +tx1^(1/par[(d)+7]) + tx2^(1/par[(d)+7]))^(2*par[(d)+7]))) +
        n*(-log(par[d+2]) - log(par[d+4])- log(par[d+6]) - 
             sum((par[1:(d)]^2)/(2*(priorscale[1:(d)]))) - 
             (priorscale[d + 1] + 1)*log(par[d+2]) - (priorscale[d + 2]/par[d+2]) -
             (priorscale[(d) + 3] + 1)*log(par[d+4]) - (priorscale[(d) + 4]/par[d+4]) -
             (priorscale[(d) + 5] + 1)*log(par[d+6]) - (priorscale[(d) + 6]/par[d+6]) +
             log(ifelse(par[d+1] > priorscale[(d) + 7] & par[d+1]< priorscale[(d) + 8], 1/(priorscale[(d) + 8]-priorscale[(d) + 7]), 0)) + 
             log(ifelse(par[d+3] > priorscale[(d) + 9] & par[d+3] < priorscale[(d) + 10], 1/(priorscale[(d) + 10] - priorscale[(d) + 9]), 0)) + 
             log(ifelse(par[d+5]> priorscale[(d) + 11] & par[d+5]< priorscale[(d) + 12], 1/(priorscale[(d) + 12] - priorscale[(d) + 11]), 0))+ 
             log(ifelse(par[(d)+7]> priorscale[(d) + 13] & par[(d)+7]< priorscale[(d) + 14], 1/(priorscale[(d) + 14] - priorscale[(d) + 13]), 0)))
    }else{
      loglik = -Inf
    }
  }
  else if(model == "HalfNormalDamagesEVD"){
    ## mu1 nonstationary
    mu1 <- par[1] + rowSums(t(t(data[,covs[[1]]])*par[2:(length(covs[[1]])+1)])) 
    mu2 <- par[length(covs[[1]]) + 2] + rowSums(t(t(data[,covs[[2]]])*par[(length(covs[[1]]) +length(covs[[2]]) + 2)]))
    mu3 <- par[(length(covs[[1]]) +length(covs[[2]]) + 3)] +  rowSums(t(t(data[,covs[[3]]])*par[(length(covs[[1]]) +length(covs[[2]]) + 4):(length(covs[[1]]) +length(covs[[2]]) + 3 + length(covs[[3]]))]))
    
    d <- (length(covs[[1]]) +length(covs[[2]]) + 3 + length(covs[[3]]))
    ## t(z1)  
    term1_tz1 <- 1 + ((par[d+1]*(data$z1 - mu1))/(par[d+2]))
    tz1_exponent <- -(1/par[d+1])
    if(par[d+1]!= 0 & all(term1_tz1 > 0)){
      tz1 <- term1_tz1^tz1_exponent
    }
    else if(par[d+1] == 0 & all(term1_tz1 > 0)){
      tz1 <- exp(-(data$z1-mu1)/par[d+2])
    }
    else{
      tz1 <- 0
    }
    ## t(x1) 
    term1_tx1 <- 1 + ((par[d+3]*(data$x1 - mu2))/(par[d+4]))
    tx1_exponent <- -(1/par[d+3])
    if(par[d+3]!= 0 & all(term1_tx1 > 0)){
      tx1 <- term1_tx1^tx1_exponent
    }
    else if(par[d+3] == 0 & all(term1_tx1 > 0)){
      tx1 <- exp(-(data$x1-mu2)/par[d+4])
    }
    else{
      tx1 <- 0
    }
    
    
    loglik_EVD <- function(y, sigma, xi){
      -log(sigma) + (xi + 1)*log(y) - y
    }
    
    if(par[d+2] > 0 & par[d+4] >0 & par[d+5] > 0 & all(tx1 != 0) & all(tz1 != 0)){
      loglik <- sum(loglik_EVD(y = tz1, sigma = par[d+2], xi = par[d+1]) + 
                      loglik_EVD(y = tx1, sigma = par[d+4], xi = par[d+3]) - 
                      log(data$x2* par[d+5] * sqrt(2*pi)) - 
                      ((log(data$x2) - mu3)^2)/(2*par[d+5]^2)) + 
        n*( - sum((par[1:(d)]^2)/(2*(priorscale[1:(d)]))) -  
              (priorscale[(d)+1] + 1)*log(par[(d)+2]) - (priorscale[(d)+2]/par[d+2]) -
              (priorscale[(d)+3] + 1)*log(par[d+4]) - (priorscale[(d)+4]/ par[d+4]) -
              (priorscale[(d)+5] + 1)*log( par[d+5]) - (priorscale[(d)+6]/ par[d+5]) +
              #log(ifelse(par[14] > priorscale[18] & par[14]< priorscale[19], 1/(priorscale[19]-priorscale[18]), 0)) +
              log(ifelse(par[d+1] > priorscale[(d)+7] & par[d+1]< priorscale[(d)+8], 1/(priorscale[(d)+8]-priorscale[(d)+7]), 0)) + 
              log(ifelse(par[d+3] > priorscale[(d)+ 9] & par[d+3] < priorscale[(d)+10], 1/(priorscale[(d)+10] - priorscale[(d)+9]), 0)))
    }else{
      loglik = -Inf
    }
  }
  return(loglik)
}
