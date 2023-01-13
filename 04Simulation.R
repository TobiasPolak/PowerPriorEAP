setwd("C:/Users/31612/Documents/R/PowerPriorPropensity")

runSimulation <- function(setting, mainFolder) {
  
  # check if the user has provided the folder to save the results 
  if(missing(mainFolder)){
    cat('You did not specify a folder, using the current folder: ', getwd(), ' to save your results.')
    mainFolder <- getwd()
  } 
  
  start.time <- Sys.time() # record the start time
  start.date <- Sys.Date() # record the start date
  
  # assign values to parameters from the input setting 
  n.sims <- setting$n.sims
  n.run <- setting$n.run
  n.settings <- setting$n.settings
  
  rct.size <- setting$rct.size
  rwd.size <- setting$rwd.size
  n.cov <- setting$n.cov
  mu.rct <- setting$mu.rct
  sigma.rct <- setting$sigma.rct
  mu.rwd <- setting$mu.rwd
  sigma.rwd <- setting$sigma.rwd
  trt.dif <- setting$trt.dif
  drift <- setting$drift
  b.coef <- setting$b.coef
  rate <- setting$rate
  mixture <- setting$mixture
  mixture.frac <- setting$mixture.frac
  wang = setting$wang
  wang2 = setting$wang2
  save_all = setting$save_all
  plotYN = setting$plotYN
  
  # if save_all is true, initialize empty lists for rates, listType, listRMSE
  if (save_all){
    rates <- list()
    listType <- list()
    listRMSE <- list()
  }
  
  # assign methods based on the value of wang2 and wang
  if (wang2) {
    methods <- c('Ignore', 'Pooling', 'ProPP', 'ProPP ATE', 'MPP', 'Wang 10%', 'Wang 20%')
  } else if (wang){
    methods <- c('Ignore', 'Pooling', 'ProPP', 'ProPP ATE', 'MPP', 'Wang 10%')
  } else{
    methods <- c('Ignore', 'Pooling', 'ProPP', 'ProPP ATE', 'MPP')
  }
  
  # initialize vectors to store results
  rate_y_rct <- rep(0, n.sims)
  rate_y_rwd <- rep(0, n.sims)
  propdelta <- rep(0, n.sims)
  propensities <- rep(0, n.sims)
  
  # Initialize matrix to store results for Type I error rate for each simulation and method
  res.type.I <- matrix(0, nrow = n.sims, ncol = length(methods))
  colnames(res.type.I) <- methods
  
  # Initialize matrix to store the mean of Type I error rate for each method for each setting
  res.type.I.mean <- matrix(0, nrow = n.settings, ncol = length(methods))
  colnames(res.type.I.mean) <- methods
  
  # Initialize matrix to store results for Root Mean Squared Error (RMSE) for each simulation and method
  RMSE <- matrix(0, nrow = n.sims, ncol = length(methods))
  colnames(RMSE) <-  methods
  
  # Initialize matrix to store the mean of RMSE for each method for each setting
  RMSE.mean <- matrix(0, nrow = n.settings, ncol = length(methods))
  colnames(RMSE.mean) <-  methods
  
  # Initialize matrix to store the mean of the posterior for each method for each simulation
  post.mean <- matrix(0, nrow = n.sims, ncol = length(methods))
  colnames(post.mean) <-  methods
  
  # Initialize matrix to store the median of the posterior for each method for each simulation
  post.median <- post.mean
  
  # Initialize matrix to store the error of the Wang estimator for each setting
  errorWang <- matrix(0, nrow = n.settings)
  
  # Initialize progress bar
  pb <- progress::progress_bar$new(
    format = " simulating [:bar] :percent eta: :eta",
    total = n.sims * n.settings,
    clear = FALSE,
    width = 60
  )       
  
  # Loop through each setting
  for (j in 1:n.settings) {
    # Create a vector of settings 
    setting.vector <- seq.int(from = setting$from,
                              to = setting$to,
                              length.out = setting$n.settings)
    # Assign the variable name
    id.vars <- 'setting.vector'
    # Update the drift or mu.rwd based on the current setting
    if (setting$variable == 'drift'){
      drift <- setting.vector[j]
    } else if(setting$variable == 'mu.rwd'){
      mu.rwd <- setting.vector[j]
    }
    # Create the data using the simulationSet() function
    data <- simulationSet(
      rct.size=rct.size*1000,
      rwd.size=rwd.size*1000,
      n.cov,
      mu.rct,
      sigma.rct,
      mu.rwd,
      sigma.rwd,
      trt.dif,
      drift,
      b.coef,
      rate,
      mixture,
      mixture.frac
    )
    # Calculate the rate of control group under H0
    rate_H0 <- mean(data$df[data$df[, "current"] == 1 & data$df[, "treat"] == 0, "prob"])
    # Loop through each simulation
    for (i in 1:n.sims) {
      # Set the seed
      set.seed(i)
      # Create the data using the simulationSet() function
      data <- simulationSet(
        rct.size,
        rwd.size,
        n.cov,
        mu.rct,
        sigma.rct,
        mu.rwd,
        sigma.rwd,
        trt.dif,
        drift,
        b.coef,
        rate,
        mixture,
        mixture.frac
      )
      
      # Estimate treatment effect using the (ProPP) method
      res <- estimateProPP(data$df, n.run = n.run, both = FALSE, rescaling = 'treated.cap')
      
      # Calculate p-value for one-sided test of treatment effect being less than rate_H0
      p.value.one.sided <- mean(res$normal.theta.sample < rate_H0)
      
      # Calculate two-sided p-value based on one-sided p-value
      p.value.two.sided <- min(2 * p.value.one.sided, 2 - 2 * p.value.one.sided)
      
      # Store the two-sided p-value, root mean squared error (RMSE), mean, and median of the estimated treatment effect
      res.type.I[i, 'Ignore'] <- p.value.two.sided
      RMSE[i, 'Ignore'] <- sqrt(mean((res$normal.theta.sample - rate_H0) ^ 2))
      post.mean[i, 'Ignore'] <- mean(res$normal.theta.sample)
      post.median[i, 'Ignore'] <- median(res$normal.theta.sample)
      
      # Repeats the process for the ProPP method
      p.value.one.sided <- mean(res$theta.post < rate_H0)
      p.value.two.sided <- min(2 * p.value.one.sided, 2 - 2 * p.value.one.sided)
      res.type.I[i, 'ProPP'] <- p.value.two.sided
      RMSE[i, 'ProPP'] <- sqrt(mean((res$theta.post - rate_H0) ^ 2))
      post.mean[i, 'ProPP'] <- mean(res$theta.post)
      post.median[i, 'ProPP'] <- median(res$theta.post)
      propensities[i] <-  res$mean.prop.rwd
      propdelta[i] <-   mean(res$propdelta)
      
      # Estimate treatment effect using the Modified Power Prior (MPP) method
      res <- estimateMPP(data$df, n.run = n.run)
      
      # Calculate p-value for one-sided test of treatment effect being less than rate
      p.value.one.sided <- mean(res$delta.MPP$p < rate_H0)
      
      # Calculate two-sided p-value based on one-sided p-value
      p.value.two.sided <- min(2 * p.value.one.sided, 2 - 2 * p.value.one.sided)
      
      # Store the two-sided p-value, root mean squared error (RMSE), mean, and median of the estimated treatment effect
      res.type.I[i, 'MPP'] <- p.value.two.sided
      RMSE[i, 'MPP'] <- sqrt(mean((res$delta.MPP$p - rate_H0) ^ 2))
      post.mean[i, 'MPP'] <- mean(res$delta.MPP$p)
      post.median[i, 'MPP'] <- median(res$delta.MPP$p)
      
      # Estimate treatment effect using the Pooling method
      res <- estimatePooled(data$df, n.run = n.run)
      
      # Calculate p-value for one-sided test of treatment effect being less than rate_H0
      p.value.one.sided <- mean(res$theta.post.pooled < rate_H0)
      
      # Calculate two-sided p-value based on one-sided p-value
      p.value.two.sided <- min(2 * p.value.one.sided, 2 - 2 * p.value.one.sided)
      
      # Store the two-sided p-value, root mean squared error (RMSE), mean, and median of the estimated treatment effect
      res.type.I[i, 'Pooling'] <- p.value.two.sided
      RMSE[i, 'Pooling'] <- sqrt(mean((res$theta.post.pooled - rate_H0) ^ 2))
      post.mean[i, 'Pooling'] <- mean(res$theta.post.pooled)
      post.median[i, 'Pooling'] <- median(res$theta.post.pooled)
      
      # Check if the "wang" variable is set to TRUE
      if (wang == TRUE) {
        # Open a file to store messages
        zz <- file("messages.Rout", open = "wt")
        # Redirect messages to the file
        sink(zz, type = "message")
        sink(zz)
        
        # Estimate propensity scores for the control group using the psrwe_est() function
        dta_ps <- psrwe_est(
          data$df[data$df[, "treat"] == 0,],
          v_covs = paste('V', 1:n.cov, sep = ""),
          v_grp = "current",
          cur_grp_level = 1,
          nstrata = 5,
          ps_method = "logistic"
        )
        
        # Use the psrwe_borrow() function to borrow information from the treatment group 
        ps_bor <- psrwe_borrow(dta_ps,
                               total_borrow = round( (setting$wang.borrowfrac  * rct.size) / 2,0),
                               method = "distance")
        
        # Use the psrwe_powerp() function to estimate treatment effect
        downloadcheck <- try({
          rst_pp <- psrwe::psrwe_powerp(
            ps_bor,
            v_outcome = "y",
            outcome_type =  "binary",
            iter = setting$wang.iter,
            warmup = setting$wang.warmup,
            chains = setting$wang.chains,
          ) 
        }, silent=TRUE)
        if (class(downloadcheck) == "try-error") {
          errorWang[j] <- errorWang[j] + 1
          next
        } else {
        }
        
        p.value.one.sided <- mean(rst_pp$Control$Overall_Samples < rate_H0)
        # calculates the one-sided p-value by taking the mean of the values in rst_pp$Control$Overall_Samples 
        # that are less than rate_H0
        
        p.value.two.sided <- min(2 * p.value.one.sided, 2 - 2 * p.value.one.sided)
        # calculates the two-sided p-value by taking the minimum of twice the one-sided p-value and 
        # twice the complement of the one-sided p-value
        
        res.type.I[i, 'Wang 10%'] <- p.value.two.sided
        # assigns the two-sided p-value to the corresponding element in the 'Wang 10%' column of the res.type.I matrix
        
        RMSE[i, 'Wang 10%'] <- sqrt(mean((rst_pp$Control$Overall_Samples - rate_H0) ^ 2))
        # calculates the root mean squared error (RMSE) between the values in rst_pp$Control$Overall_Samples and rate_H0,
        # and assigns the value to the corresponding element in the 'Wang 10%' column of the RMSE matrix
        
        post.mean[i, 'Wang 10%'] <- mean(rst_pp$Control$Overall_Samples)
        # assigns the mean of the values in rst_pp$Control$Overall_Samples to the corresponding element in the 'Wang 10%' column of the post.mean matrix
        
        post.median[i, 'Wang 10%'] <- median(rst_pp$Control$Overall_Samples)
        # assigns the median of the values in rst_pp$Control$Overall_Samples to the corresponding element in the 'Wang 10%' column of the post.median matrix
        
        # if wang2 is TRUE, the code inside the if statement will be executed. This block of code is similar to the previous one, but with the difference of using 'Wang 20%' instead of 'Wang 10%' and using ps_bor object to perform the power calculation
        
        if (wang2 == TRUE){
          ps_bor <- psrwe_borrow(dta_ps,
                                 total_borrow = round( (0.2  * rct.size) / 2,0),
                                 method = "distance")
          
          rst_pp <- psrwe::psrwe_powerp(
            ps_bor,
            v_outcome = "y",
            outcome_type = "binary",
            iter = setting$wang.iter,
            warmup = setting$wang.warmup,
            chains = setting$wang.chains,
          )
          
          # rate_H0 <-
          #   mean(data$df[data$df[, "current"] == 1 &
          #                  data$df[, "treat"] == 0, "prob"])
          p.value.one.sided <-
            mean(rst_pp$Control$Overall_Samples < rate_H0)
          p.value.two.sided <-
            min(2 * p.value.one.sided, 2 - 2 * p.value.one.sided)
          res.type.I[i, 'Wang 20%'] <- p.value.two.sided
          RMSE[i, 'Wang 20%'] <-
            sqrt(mean((
              rst_pp$Control$Overall_Samples - rate_H0
            ) ^ 2))
          post.mean[i, 'Wang 20%'] <- mean(rst_pp$Control$Overall_Samples)
          post.median[i, 'Wang 20%'] <- median(rst_pp$Control$Overall_Samples)
        }
        
        sink()
        # closes the connection to the sink, which is used to redirect output to a file or other connection
        
        closeAllConnections()
        # closes all open connections
      }
      
      rate_y_rct[i] <- mean(data$df[data$df[, "current"] == 1 & data$df[, "treat"] == 0, "y"])
      # calculates mean of the 'y' column where "current" column is equal to 1 and "treat" column is equal to 0 and assigns the value to rate_y_rct[i]
      
      rate_y_rwd[i] <- mean(data$df[data$df[, "current"] == 0 & data$df[, "treat"] == 0, "y"])
      # calculates mean of the 'y' column where "current" column is equal to 0 and "treat" column is equal to 0 and assigns the value to rate_y_rwd[i]
      
      pb$tick()
      # increments the progress bar
    }
    res.type.I.mean[j,] <- colMeans(res.type.I < 0.05, na.rm=T)
    # calculates the mean of columns of res.type.I matrix where elements are less than 0.05 and assigns the value to the corresponding element in the jth row of res.type.I.mean matrix. The missing values are removed before calculating the mean
    
    RMSE.mean[j,] <- colMeans(RMSE, na.rm=T)
    # calculates the mean of all columns of RMSE matrix and assigns the value to the corresponding element in the jth row of RMSE.mean matrix. The missing values are removed before calculating the mean
    
    if (save_all){
      listRMSE <- list.append(listRMSE, list(RMSE = RMSE, iter = j, value = setting.vector[j]))
      listType <- list.append(listType, list(TypeI = res.type.I, iter = j, value = setting.vector[j]))
      rates <- list.append(rates, list( iter = j, 
                                        value = setting.vector[j],
                                        rate_y_rct = rate_y_rct, 
                                        rate_y_rwd = rate_y_rwd,
                                        post.mean = post.mean,
                                        post.median = post.median,
                                        propensities = propensities,
                                        propdelta = propdelta
      ))
    } else{}
    
    write.csv2( 
      res.type.I.mean,
      paste0(start.date, '_res.type.I.mean_', setting$title, j, '.csv')
    )
    
    write.csv2( 
      RMSE.mean,
      paste0(start.date, '_RMSE.mean_', setting$title, j, '.csv')
    )
  }
  res.type.I.mean <- cbind(setting.vector, res.type.I.mean)
  RMSE.mean <- cbind(setting.vector, RMSE.mean)
  
  mainDir <- mainFolder # Assign mainFolder to mainDir
  mainDir <- "C:/Users/31612/Documents/R/PowerPriorPropensity/Simulations" # Re-assign mainDir to specific file path
  folder_name <- paste0('Simulation_', start.date) # Create folder_name variable with "Simulation_" and start.date
  subDir <- paste0("/", folder_name) # Create subDir variable with subfolder name
  dir.create(file.path(mainDir, subDir), showWarnings = FALSE) # Create new directory
  subfolder <- paste0('/',setting$title) # Create subfolder variable with setting$title
  dir.create(file.path(mainDir, subDir, subfolder), showWarnings = FALSE) # Create subfolder within directory
  
  setwd(file.path(file.path(mainDir, subDir, subfolder))) # Set working directory to subfolder
  
  # Plot and save Type I error and RMSE CSV
  
  write.csv2( 
    res.type.I.mean,
    paste0(start.date, '_res.type.I.mean_', setting$title, '.csv')
  )
  
  write.csv2( 
    RMSE.mean,
    paste0(start.date, '_RMSE.mean_', setting$title, '.csv')
  )
  
  # Plot and save Type I error and RMSE graphs
  
  if (plotYN){
    plot.df.type.I <-
      melt(
        as.data.frame(res.type.I.mean),
        id.vars = id.vars,
        variable.name = 'Method'
      )
    
    ggplot(data = plot.df.type.I,
           aes(
             x = setting.vector,
             y = value,
             group = Method,
             color = Method
           )) +
      geom_line() + theme_bw() + geom_point() + labs(title = paste0(setting$title, " Simulation Type I error"),
                                                     x = setting$variable, y = "Type I error rate")
    ggsave(filename = paste0(setting$title, " Simulation Type I error", ".eps"), device = 'eps', dpi = 600, width = 30, height = 20, units = "cm")
    ggsave(filename = paste0(setting$title, " Simulation Type I error", ".png"), device = 'png', dpi = 600, width = 30, height = 20, units = "cm")
    
    
    plot.df.RMSE <-
      melt(as.data.frame(RMSE.mean),
           id.vars = id.vars,
           variable.name = "Method")
    ggplot(data = plot.df.RMSE,
           aes(
             x = setting.vector,
             y = value,
             group = Method,
             color = Method
           )) +
      
      geom_line() + theme_bw() + geom_point() + labs(title = paste0(setting$title," Simulation RMSE"),
                                                     x = setting$variable, y = "RMSE")
    ggsave(filename = paste0(setting$title, " Simulation RMSE", ".eps"), device = 'eps', dpi = 600, width = 30, height = 20, units = "cm")
    ggsave(filename = paste0(setting$title, " Simulation RMSE", ".png"), device = 'png', dpi = 600, width = 30, height = 20, units = "cm")
  } else{}
  
  # Record end time, calculate time spent, and save data
  end.time <- Sys.time()  # Record the current time and assign it to end.time variable
  end.date <- Sys.Date()  # Record the current date and assign it to end.date variable
  timespent <- end.time - start.time  # Calculate the time spent by subtracting start.time from end.time
  print(timespent)  # Print the time spent
  if (save_all){
    saveRDS(listRMSE, file=paste0(start.date,'_', setting$title, '_' , setting$variable ,"listRMSE.RData"))  # Save listRMSE variable in .RData format
    saveRDS(listType, file=paste0(start.date,'_', setting$title, '_' , setting$variable ,"listRMSE.RData"))  # Save listType variable in .RData format
    saveRDS(rates, file=paste0(start.date,'_', setting$title, '_' , setting$variable ,"rates.RData"))  # Save rates variable in .RData format
    save.image(file = paste0(start.date, setting$variable ,".RData"))  # Save the current R session in .RData format
  }
  information_path <- paste0(start.date,'_', setting$title, '_', 'information.txt')  # Create a variable for the information file path
  write.table(setting, file = information_path)  # Write setting variable to the information file
  write.table(paste('Start time: ', start.time, 'End time: ', end.time, 'Time spent: ', timespent), file = information_path, append=TRUE, col.names = FALSE)  # Append the start time, end time and time spent to the information file
  write.table(paste('Error Wang: ', errorWang/n.sims), file = information_path, append=TRUE, col.names = FALSE)  # Append the Error Wang value to the information file
}
###########################
###########################
###########################





