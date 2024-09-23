### Function for Posterior in a hierarchical GEV/trivariate/GEV with log-Normal damages model with  damages, windspeed and pressure as variables in the hierarchy.
## Input: # priorscale: hyperparameter choices,
          # data: to be inputted
          # covs: covariates common to each layer of hierarchy
          # model: either "Hier_GEV" or "Tri_GEV" or "HalfNormalDamagesEVD" (depending on the model being used)
## Output: # unnormalized posterior from the log-likelihood and priors

logpost_factory <- function(priorscale, data, covs) function(par,  model){
  if(model == "Hier_GEV"){
  n <- nrow(data)  
  covs_ind <- which(colnames(data) %in% covs)
  d <- length(covs)
  ## nonstationary location parameters
  mu1 <- par[1] + rowSums(t(t(data[,covs_ind[1:(d-2)]])*par[2:(d-1)])) 
  mu2 <- par[d] + rowSums(t(t(data[,covs_ind[c(d,1:(d-2))]])*par[(d+1):(2*d-1)]))
  mu3 <- par[2*d] +  rowSums(t(t(data[,covs_ind[c((d-1):d,1:(d-2))]])*par[(2*d + 1):(3*d)]))
  
  ## t(z1)  
  term1_tz1 <- 1 + ((par[d*3+1]*(data$z1 - mu1))/(par[d*3+2]))
  tz1_exponent <- -(1/par[d*3+1])
  if(par[d*3+1]!= 0 & all(term1_tz1 > 0)){
    tz1 <- term1_tz1^tz1_exponent
  }
  else if(par[d*3+1] == 0 & all(term1_tz1 > 0)){
    tz1 <- exp(-(data$z1-mu1)/par[d*3+2])
  }
  else{
    tz1 <- 0
  }
  ## t(x1) 
  term1_tx1 <- 1 + ((par[d*3+3]*(data$x1 - mu2))/(par[d*3+4]))
  tx1_exponent <- -(1/par[d*3+3])
  if(par[d*3+3]!= 0 & all(term1_tx1 > 0)){
    tx1 <- term1_tx1^tx1_exponent
  }
  else if(par[d*3+3] == 0 & all(term1_tx1 > 0)){
    tx1 <- exp(-(data$x1-mu2)/par[d*3+4])
  }
  else{
    tx1 <- 0
  }
  ## t(x2) 
  term1_tx2 <- 1 + ((par[d*3+5]*(data$x2 - mu3))/(par[d*3+6]))
  tx2_exponent <- -(1/par[d*3+5])
  if(par[d*3+5]!= 0 & all(term1_tx2 > 0)){
    tx2 <- term1_tx2^tx2_exponent
  }
  else if(par[d*3+5] == 0 & all(term1_tx2 > 0)){
    tx2 <- exp(-(data$x2-mu3)/par[d*3+6])
  }
  else{
    tx2 <- 0
  }
  loglik_EVD <- function(y, sigma, xi){
    -log(sigma) + (xi + 1)*log(y) - y
  }
  
  if(par[d*3+2] > 0 & par[d*3+4] >0 & par[d*3+6] > 0 & all(tz1 != 0) & all(tx1 != 0) & all(tx2 != 0)){
    loglik <- sum(loglik_EVD(y = tz1, sigma = par[d*3+2], xi = par[d*3+1]) + 
                    loglik_EVD(y = tx1, sigma = par[d*3+4], xi = par[d*3+3]) + 
                    loglik_EVD(y = tx2, sigma = par[d*3+6], xi = par[d*3+5])) +
      n*( - sum((par[1:(d*3)]^2)/(2*(priorscale[1:(d*3)])))- 
            (priorscale[d*3 + 1] + 1)*log(par[d*3+2]) - (priorscale[d*3 + 2]/par[d*3+2]) -
            (priorscale[(d*3) + 3] + 1)*log(par[d*3+4]) - (priorscale[(d*3) + 4]/par[d*3+4]) -
            (priorscale[(d*3) + 5] + 1)*log(par[d*3+6]) - (priorscale[(d*3) + 6]/par[d*3+6]) +
            log(ifelse(par[d*3+1] > priorscale[(d*3) + 7] & par[d*3+1]< priorscale[(d*3) + 8], 1/(priorscale[(d*3) + 8]-priorscale[(d*3) + 7]), 0)) + 
            log(ifelse(par[d*3+3] > priorscale[(d*3) + 9] & par[d*3+3] < priorscale[(d*3) + 10], 1/(priorscale[(d*3) + 10] - priorscale[(d*3) + 9]), 0)) + 
            log(ifelse(par[d*3+5]> priorscale[(d*3) + 11] & par[d*3+5]< priorscale[(d*3) + 12], 1/(priorscale[(d*3) + 12] - priorscale[(d*3) + 11]), 0)))
    
  }else{
    loglik = -Inf
  }}
  else if(model == "Tri_GEV"){
    n <- nrow(data)  
    covs_ind <- which(colnames(data) %in% covs)
    d <- length(covs) 
    ## mu1 nonstationary
    mu1 <- par[1] + rowSums(t(t(data[,covs_ind[1:(d-2)]])*par[2:(d-1)])) 
    mu2 <- par[d] + rowSums(t(t(data[,covs_ind[c(d,1:(d-2))]])*par[(d+1):(2*d-1)]))
    mu3 <- par[2*d] +  rowSums(t(t(data[,covs_ind[c((d-1):d,1:(d-2))]])*par[(2*d + 1):(3*d)]))
    
    ## t(z1)  
    term1_tz1 <- 1 + ((par[d*3+1]*(data$z1 - mu1))/(par[d*3+2]))
    tz1_exponent <- -(1/par[d*3+1])
    if(par[d*3+1]!= 0 & all(term1_tz1 > 0)){
      tz1 <- term1_tz1^tz1_exponent
    }
    else if(par[d*3+1] == 0 & all(term1_tz1 > 0)){
      tz1 <- exp(-(data$z1-mu1)/par[d*3+2])
    }
    else{
      tz1 <- 0
    }
    ## t(x1) 
    term1_tx1 <- 1 + ((par[d*3+3]*(data$x1 - mu2))/(par[d*3+4]))
    tx1_exponent <- -(1/par[d*3+3])
    if(par[d*3+3]!= 0 & all(term1_tx1 > 0)){
      tx1 <- term1_tx1^tx1_exponent
    }
    else if(par[d*3+3] == 0 & all(term1_tx1 > 0)){
      tx1 <- exp(-(data$x1-mu2)/par[d*3+4])
    }
    else{
      tx1 <- 0
    }
    ## t(x2) 
    term1_tx2 <- 1 + ((par[d*3+5]*(data$x2 - mu3))/(par[d*3+6]))
    tx2_exponent <- -(1/par[d*3+5])
    if(par[d*3+5]!= 0 & all(term1_tx2 > 0)){
      tx2 <- term1_tx2^tx2_exponent
    }
    else if(par[d*3+5] == 0 & all(term1_tx2 > 0)){
      tx2 <- exp(-(data$x2-mu3)/par[d*3+6])
    }
    else{
      tx2 <- 0
    }
    loglik_EVD <- function(y, sigma, xi){
      -log(sigma) + (xi + 1)*log(y) - y
    }
    if(par[d*3+2] > 0 & par[d*3+4] >0 & par[d*3+6] > 0 & all(tz1 != 0) & all(tx1 != 0) & all(tx2 != 0) & par[d*3+7] > 0.05 & par[d*3+7] < 1){
      loglik <- sum((par[d*3+1]+1)*log(tz1) + (par[d*3+3]+1)*log(tx1) + (par[d*3+5]+1)*log(tx2) +
                      ((1/par[d*3+7]) - 1)*log(tz1) + ((1/par[d*3+7])-1)*log(tx1) + ((1/par[d*3+7])-1)*log(tx2) - 
                      (tx1^(1/par[d*3+7]) + tx2^(1/par[d*3+7]) + tz1^(1/par[d*3+7]))^par[d*3+7] + 
                      (par[(d*3)+7]-3)*log(tz1^(1/par[(d*3)+7]) + tx1^(1/par[(d*3)+7]) + tx2^(1/par[(d*3)+7])) + 
                      log((((1-par[(d*3)+7])*(2-par[(d*3)+7]))/par[(d*3)+7]) - ((1-par[(d*3)+7])/par[(d*3)+7])*(tz1^(1/par[(d*3)+7]) +tx1^(1/par[(d*3)+7]) + tx2^(1/par[(d*3)+7]))^par[(d*3)+7] + 
                            (tz1^(1/par[(d*3)+7]) +tx1^(1/par[(d*3)+7]) + tx2^(1/par[(d*3)+7]))^(2*par[(d*3)+7]))) +
        n*(-log(par[d*3+2]) - log(par[d*3+4])- log(par[d*3+6]) - 
             sum((par[1:(d*3)]^2)/(2*(priorscale[1:(d*3)]))) - 
             (priorscale[d*3 + 1] + 1)*log(par[d*3+2]) - (priorscale[d*3 + 2]/par[d*3+2]) -
             (priorscale[(d*3) + 3] + 1)*log(par[d*3+4]) - (priorscale[(d*3) + 4]/par[d*3+4]) -
             (priorscale[(d*3) + 5] + 1)*log(par[d*3+6]) - (priorscale[(d*3) + 6]/par[d*3+6]) +
             log(ifelse(par[d*3+1] > priorscale[(d*3) + 7] & par[d*3+1]< priorscale[(d*3) + 8], 1/(priorscale[(d*3) + 8]-priorscale[(d*3) + 7]), 0)) + 
             log(ifelse(par[d*3+3] > priorscale[(d*3) + 9] & par[d*3+3] < priorscale[(d*3) + 10], 1/(priorscale[(d*3) + 10] - priorscale[(d*3) + 9]), 0)) + 
             log(ifelse(par[d*3+5]> priorscale[(d*3) + 11] & par[d*3+5]< priorscale[(d*3) + 12], 1/(priorscale[(d*3) + 12] - priorscale[(d*3) + 11]), 0))+ 
             log(ifelse(par[(d*3)+7]> priorscale[(d*3) + 13] & par[(d*3)+7]< priorscale[(d*3) + 14], 1/(priorscale[(d*3) + 14] - priorscale[(d*3) + 13]), 0)))
    }else{
      loglik = -Inf
    }}
  else if(model == "HalfNormalDamagesEVD"){
    n <- nrow(data)  
    covs_ind <- which(colnames(data) %in% covs)
    d <- length(covs) 
    
    ## mu1 nonstationary
    mu1 <- par[1] + rowSums(t(t(data[,covs_ind[1:(d-2)]])*par[2:(d-1)])) 
    mu2 <- par[d] + rowSums(t(t(data[,covs_ind[c(d,1:(d-2))]])*par[(d+1):(2*d-1)]))
    mu3 <- par[2*d] +  rowSums(t(t(data[,covs_ind[c((d-1):d,1:(d-2))]])*par[c((2*d+2),(2*d + 1),(2*d+3):(3*d))]))
    
    ## t(z1)  
    term1_tz1 <- 1 + ((par[d*3+1]*(data$z1 - mu1))/(par[d*3+2]))
    tz1_exponent <- -(1/par[d*3+1])
    if(par[d*3+1]!= 0 & all(term1_tz1 > 0)){
      tz1 <- term1_tz1^tz1_exponent
    }
    else if(par[d*3+1] == 0 & all(term1_tz1 > 0)){
      tz1 <- exp(-(data$z1-mu1)/par[d*3+2])
    }
    else{
      tz1 <- 0
    }
    ## t(x1) 
    term1_tx1 <- 1 + ((par[d*3+3]*(data$x1 - mu2))/(par[d*3+4]))
    tx1_exponent <- -(1/par[d*3+3])
    if(par[d*3+3]!= 0 & all(term1_tx1 > 0)){
      tx1 <- term1_tx1^tx1_exponent
    }
    else if(par[d*3+3] == 0 & all(term1_tx1 > 0)){
      tx1 <- exp(-(data$x1-mu2)/par[d*3+4])
    }
    else{
      tx1 <- 0
    }
    
    
    loglik_EVD <- function(y, sigma, xi){
      -log(sigma) + (xi + 1)*log(y) - y
    }
    
    if(par[d*3+2] > 0 & par[d*3+4] >0 & par[d*3+5] > 0 & all(tx1 != 0) & all(tz1 != 0)){
      loglik <- sum(loglik_EVD(y = tz1, sigma = par[d*3+2], xi = par[d*3+1]) + 
                      loglik_EVD(y = tx1, sigma = par[d*3+4], xi = par[d*3+3]) - 
                      log(data$x2* par[d*3+5] * sqrt(2*pi)) - 
                      ((log(data$x2) - mu3)^2)/(2*par[d*3+5]^2)) + 
        n*( - sum((par[1:(d*3)]^2)/(2*(priorscale[1:(d*3)]))) -  
              (priorscale[(d*3)+1] + 1)*log(par[(d*3)+2]) - (priorscale[(d*3)+2]/par[d*3+2]) -
              (priorscale[(d*3)+3] + 1)*log(par[d*3+4]) - (priorscale[(d*3)+4]/ par[d*3+4]) -
              (priorscale[(d*3)+5] + 1)*log( par[d*3+5]) - (priorscale[(d*3)+6]/ par[d*3+5]) +
              #log(ifelse(par[14] > priorscale[18] & par[14]< priorscale[19], 1/(priorscale[19]-priorscale[18]), 0)) +
              log(ifelse(par[d*3+1] > priorscale[(d*3)+7] & par[d*3+1]< priorscale[(d*3)+8], 1/(priorscale[(d*3)+8]-priorscale[(d*3)+7]), 0)) + 
              log(ifelse(par[d*3+3] > priorscale[(d*3)+ 9] & par[d*3+3] < priorscale[(d*3)+10], 1/(priorscale[(d*3)+10] - priorscale[(d*3)+9]), 0)))
    }else{
      loglik = -Inf
    }
  }
  return(loglik)
}