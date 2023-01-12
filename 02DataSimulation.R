## Input parameters
simulationSet <- function (rct.size, #Sample size of the RCT
                           rwd.size, # Sample size of the RWD
                           n.cov, # Number of covariates
                           mu.rct, # Mean vector of the covariates in the RCT
                           sigma.rct, # Variance of the RCT
                           mu.rwd, #Mean vector of the covariates in the RWD
                           sigma.rwd, # Variance of RWD
                           trt.dif, # Treatment effect
                           drift, # Drift between RCT and RWD
                           b.coef, # Coefficients of the covariates
                           rate, # Rate of missing data
                           mixture, # Whether to use mixture distribution
                           mixture.frac) {
  
  # Total sample size
  n <- rct.size + rwd.size
  
  # Repeat the mean and standard deviation of the covariates for the number of covariates
  mu.rct <- rep(mu.rct, n.cov)
  sigma.rct <- rep(sigma.rct, n.cov)
  mu.rwd <- rep(mu.rwd, n.cov)
  sigma.rwd <- rep(sigma.rwd, n.cov)
  
  # Repeat the coefficient of the covariates for the number of covariates
  if (length(b.coef) != n.cov){
    b.coef <- rep(b.coef[1], n.cov)
  }
  
  ## Simulate RCT and RWD
  current <- rep(1, n) # Initialize all subjects as current
  rwd <- sort(sample(1:n, rwd.size, replace = F))
  rct <- c(1:n)[-rwd]
  current[rwd] <- 0     #RCT = 1, RWD = 0
  noncurrent <- 1-current
  treat <- rep(1, n)    #TREATMENT = 1, CONTROL = 0
  trt.rwd <- sort(sample(rwd, round(length(rwd) / 2), replace = F))
  treat[trt.rwd] <- 0       #control group of RWD
  trt.rct <- sort(sample(rct, round(length(rct) / 2), replace = F))
  treat[trt.rct] <- 0  #control group of RCT
  
  ## Simulate the covariates
  V <- matrix(NA, nrow = n, ncol = n.cov)
  V <- matrix(rnorm(n.cov * n, mean = mu.rct, sd = sigma.rct), n, n.cov, byrow = T)
  
  ## Generate perturbation via mixture distribution
  if (mixture == TRUE) {
    if (is.na(mixture.frac)) {
      mixture.frac = 2
    }
    mixt.rwd <- sample(rwd, round(length(rwd)) / mixture.frac, replace = F)
    ## Randomly select a fraction of RWD subjects to use the mean and variance of the RWD distribution
    V[mixt.rwd, ] <- matrix(rnorm(n.cov * length(mixt.rwd), mean = mu.rwd, sd = sigma.rwd),
                            length(mixt.rwd), n.cov, byrow = T)
  }
  
  colnames(V) <- paste0('V', 1:n.cov)
  
  ## Generate outcome variable using a logistic model with treatment effect, drift, and covariates
  if(missing(rate)){
    rate <- 0 #rate is zero if not explicitly defined.
  } else {}
  
  z <- t(rate + trt.dif * treat + drift * noncurrent +  b.coef %*% t(V))
  y <- rbinom(n = n, 1, p = plogis(z)) 
  
  mean_z_rwd <- rate + trt.dif * treat + drift * noncurrent +  sum(b.coef*mu.rwd)
  mean_z_rct <- rate + trt.dif * treat + drift * noncurrent +  sum(b.coef*mu.rct)
  
  var_z_rwd <- sum((b.coef^2)*sigma.rwd^2)
  var_z_rct <- sum((b.coef^2)*sigma.rct^2)
  
  prob <- plogis(z)
  colnames(prob) <- 'prob'  
  colnames(z) <- 'Z'
  
  ## Create resulting data frame
  df <- as.data.frame(cbind(y, treat, current, V,z = z,prob=prob))
  df.h <- df[which(df$treat == 0), ]
  rct <- which(df.h$current == 1) 
  rwd <- which(df.h$current == 0)
  sum.rct <- sum(df.h$y[rct])
  sum.rwd <- sum(df.h$y[rwd])
  rate.rct <- sum(df.h$y[rct]) / length(rct) 
  rate.rwd <- sum(df.h$y[rwd]) / length(rwd)
  
  info <-
    c(
      'N_RCT' = length(rct),
      'N_RWD' = length(rwd),
      'Event_RCT' = sum.rct,
      'Event_RWD' = sum.rwd,
      'Rate_RCT' = round(rate.rct, 3),
      'Rate_RWD' =  round(rate.rwd, 3),
      'Cov_RCT' = mu.rct,
      'Cov_RWD' = mu.rwd,
      'trt.dif' = trt.dif,
      'drift' = drift,
      'beta.coef' = b.coef
    )
  result <- list(df = df, info = info)
  return(result)
}



                            