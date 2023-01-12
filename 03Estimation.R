
# This function estimates the 'ProPP' method which is propensity-score weighted power-prior in binomial model.
# It takes 4 inputs as arguments:
#   1. dataframe 'df' from 02 DataSimulation,
#   2. number of runs 'n.run',
#   3. binary 'both' for including propensities for RWD AND RCT,
#   4. input 'rescaling' for ways to rescale the propensities.
#     estimateMPP:
#   From the NPP package: estimates binominal model with modified power prior.

estimateProPP <- function (df, n.run, both, rescaling){
  
  #select only control treatments
  df <- df[which(df$treat== 0),]
  
  #count number of covariates
  n.cov <- sum(grepl("V", colnames(df)))
  
  #define a formula to be used in glm function
  covariates <- paste0("V", 1:n.cov)
  predformula <- as.formula(paste("current ~ ", paste(covariates, collapse='+')))
  
  #fit the glm model
  pred <- glm(predformula, family = binomial, data=df)
  
  #get the propensities from the fitted model
  propensity <- pred$fitted.values
  
  #add the propensities to the dataframe
  df <- cbind(df, propensity)
  
  #get the indices of RCT and RWD
  rct <- which(df$current==1)
  rwd <- which(df$current==0)
  
  ### For the calculation of delta
  ### Perform 
  
  # Check if n.run is missing, if so use default value of 10,000
  if(missing(n.run)){
    n.run <- 10.000
    print('Using default value n.run of 10.000')
  } else {}
  
  # Create delta variable with random values between 0 and 1
  delta = runif(n.run, 0, 1)
  
  ## Rescaling allows to rescale the propensity score weights.
  # Check if rescaling is missing, if so do not apply any rescaling
  if (missing(rescaling)){
    print('Rescaling is missing. No rescaling applied')
    rescaling <- 'N/A'
  } else {
    # print(rescaling)
  }
  
  # Check the value of rescaling and apply corresponding rescaling method
  if (rescaling == 'total'){ #rescale to the total number of patients in the RWD
    df$propensity[rwd] <- df$propensity[rwd] * length(rwd) / sum.prop.rwd
  } else if (rescaling == 'median'){ #rescale relative to the median propensity in the RCT
    df$propensity[rwd] <- df$propensity[rwd] / median(df$propensity[rct])
  } else if (rescaling == 'mean'){ #rescale relative to the mean propensity in the RCT
    df$propensity[rwd] <- df$propensity[rwd] / mean(df$propensity[rct])
  } else if (rescaling == 'max'){ #rescale relative to the largest propensity in the RCT
    df$propensity[rwd] <- df$propensity[rwd] / max(df$propensity[rct])
  } else if (rescaling == 'treated'){ # ATT
    df$propensity[rwd] <- df$propensity[rwd] / (1 - df$propensity[rwd])
  } else if (rescaling == 'treated.cap'){ #rescale relative to mean, maximize at 1
    df$propensity[rwd] <- pmin(1, df$propensity[rwd] / (1 - df$propensity[rwd]))
  }else{}
  
  ## Both allows (if TRUE) to weight both RCT and RWD patients with their propensity.
  if (missing(both)){
    both <- FALSE
    print('Using propensities of RWD only')
  } else {}
  if (both == FALSE){
    df$propensity[rct] <- 1
  } else {}
  #sum of RCT, RWD, propensities of RCT, propensities of RWD, mean of propensities of RCT, mean of propensities of RWD
  sum.rct <- sum(df$y[rct])
  sum.rwd <- sum(df$y[rwd])
  sum.prop.rct <- sum(df$propensity[rct])
  sum.prop.rwd <- sum(df$propensity[rwd])
  mean.prop.rct <- mean(df$propensity[rct])
  mean.prop.rwd <- mean(df$propensity[rwd])
  
  ## Acceptance sampling
  # perform acceptance sampling for both = TRUE, and for both = FALSE
  if(both == TRUE){
    shape1.numer <- delta*sum(df$propensity[rwd] * df$y[rwd]) + sum(df$propensity[rct] * df$y[rct]) + 1
    shape2.numer <- delta*sum(df$propensity[rwd] * (1-df$y[rwd])) + sum(df$propensity[rct] * (1-df$y[rct])) + 1
    shape1.denom <- delta*sum(df$propensity[rwd] * df$y[rwd]) + 1
    shape2.denom <- delta*sum(df$propensity[rwd] *(1-df$y[rwd]))+ 1 
  }  else{
    shape1.numer <- delta*sum(df$propensity[rwd] * df$y[rwd]) + sum(df$y[rct]) + 1
    shape2.numer <- delta*sum(df$propensity[rwd] * (1-df$y[rwd])) + sum(1-df$y[rct]) + 1
    shape1.denom <- delta*sum(df$propensity[rwd] * df$y[rwd]) + 1
    shape2.denom <- delta*sum(df$propensity[rwd] * (1-df$y[rwd])) + 1
  }
  
  lnumer <- lbeta(shape1.numer, shape2.numer)
  ldenom <- lbeta(shape1.denom, shape2.denom)
  scale.post <- exp(lnumer - ldenom - max(lnumer-ldenom))
  
  u = runif(n.run, 0 ,1)
  accept <- u < scale.post
  n.delta <- sum(accept)
  
  # for subsequent conditional of theta, only select delta[accept]
  lnumer <- lnumer[accept]
  ldenom <- ldenom[accept]
  scale.post <- exp(lnumer - ldenom - max(lnumer-ldenom))
  
  prop.delta.sample <- sample(df$propensity[rwd], n.run, replace = T) * sample(delta[accept], n.run, replace = T) # om delta * propensity te samplen
  prop.theta.sample <- rbeta(n.delta, shape1 = shape1.numer[accept], shape2 = shape2.numer[accept])
  normal.theta.sample <- rbeta(n.delta, shape1 = sum(df$y[rct]) + 1, shape2 = sum(1 - df$y[rct]) + 1) 
  prop.theta.sample <- rbeta(n.delta, shape1 = shape1.numer[accept], shape2 = shape2.numer[accept])
  normal.theta.sample <- rbeta(n.delta, shape1 = sum(df$y[rct]) + 1, shape2 = sum(1 - df$y[rct]) + 1) 

  result <- list(delta.post = delta[accept], 
                 rct = rct,
                 rwd = rwd,
                 sum.rct = sum.rct, 
                 total.rct = length(rct),  
                 sum.rwd = sum.rwd, 
                 total.rwd = length(rwd), 
                 ProPP.df = df, 
                 ESS = shape1.numer + shape2.numer,
                 mean.prop.rct = mean(df$propensity[rct]), 
                 mean.prop.rwd = mean(df$propensity[rwd]), 
                 propdelta = prop.delta.sample, 
                 theta.post = prop.theta.sample,
                 normal.theta.sample = normal.theta.sample,
                 both = both,
                 orig.prop = pred$fitted.values,
                 rescaling = rescaling) 
  
  return(result)
  
}

estimateMPP <- function (df, n.run){
  # remove all observations from the dataframe that are not control treatments 
  df <- df[which(df$treat== 0),]
  
  # check if n.run is missing and set default value if it is 
  if(missing(n.run)){
    n.run <- 10.000
    print('Using default value n.run of 10.000')
  } else {}
  
  # get the indices of RCT and RWD
  rct <- which(df$current==1)
  rwd <- which(df$current==0)
  sum.rct <- sum(df$y[rct])
  sum.rwd <- sum(df$y[rwd])
  
  # estimate delta using BerNPP_MCMC function from NPP package
  delta.MPP <- NPP::BerNPP_MCMC(Data.Hist = c(length(rwd), sum.rwd), Data.Cur = c(length(rct), sum.rct),
                                prior = list(p.alpha = 1, p.beta = 1, delta.alpha = 1, delta.beta = 1),
                                MCMCmethod = 'RW', rw.logit.delta = 1, nsample = n.run,
                                control.mcmc = list(delta.ini = NULL, burnin = 2000, thin = 5))
  
  # store the results in a list
  result <- list(delta.MPP = delta.MPP)
  return(result)
}

############################################################


estimatePooled <- function (df, n.run){
  
  #select only control treatments
  df <- df[which(df$treat== 0),]
  #get the indices of RCT and RWD
  rct <- which(df$current==1)
  rwd <- which(df$current==0)
  
  # check if n.run is missing and set default value if it is 
  if(missing(n.run)){
    n.run <- 10.000
    print('Using default value n.run of 10.000')
  } 
  
  # create variables for the shape parameters of the beta distribution
  shape1 <- sum(df$y[rwd]) + sum(df$y[rct]) + 1
  shape2 <- sum((1-df$y[rwd])) + sum(1-df$y[rct]) + 1
  
  # generate samples from the beta distribution
  prop.theta.sample <- rbeta(n.run, shape1 = shape1, shape2 = shape2)
  
  # store the results in a list
  result <- list(theta.post.pooled = prop.theta.sample) 
  
  return(result)
}
