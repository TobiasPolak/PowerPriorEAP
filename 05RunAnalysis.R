# Clear the environment, close all connections, set working directory and source necessary files
rm(list=ls()) # Clear all objects in the environment
closeAllConnections() # Close all open connections
setwd("C:/Users/31612/Documents/R/PowerPriorPropensity") # Set the working directory
source("01LoadInstall.R") # Source the Load and Installation R Script
source("02DataSimulation.R") # Source the Data Simulation R Script
source("03Estimation.R") # Source the Estimation R Script
source("04Simulation.R") # Source the Simulation R Script

mainFolder <- "C:/Users/31612/Documents/R/PowerPriorPropensity/Simulations" # Set main folder for simulation results

# Define several settings for the simulation
setting1 <- list(
  n.settings = 7,
  rct.size = 800,
  rwd.size = 800,
  n.cov = 5,
  mu.rct = 0,
  sigma.rct = 1,
  mu.rwd = 0,
  sigma.rwd = 1,
  trt.dif = 0,
  drift = 0,
  b.coef = 0.2,
  rate = 0,
  mixture = TRUE,
  mixture.frac = 2,
  wang = FALSE,
  wang2 = FALSE,
  wang.iter = 2500,
  wang.warmup = 1000,
  wang.chains = 3,
  wang.borrowfrac = 0.1,
  n.sims = 12500,
  n.run =  12500,
  title = 'Setting 1',
  variable = 'drift',
  from = -0.5, 
  to = 0,
  save_all = TRUE,
  plotYN = FALSE)

setting11 <- setting1
setting11[['rwd.size']] <- 5*setting11[['rct.size']]
setting11[['title']] <- 'Setting 11'

setting12 <- setting1
setting12[['rct.size']] <- 2*setting11[['rct.size']]
setting12[['title']] <- 'Setting 12'

setting13 <- setting1
setting13[['n.cov']] <- 10
setting13[['title']] <- 'Setting 13'

setting2 <- setting1
setting2[['title']] <- 'Setting 2'
setting2[['variable']] <- 'mu.rwd'

setting21 <- setting2
setting21[['rwd.size']] <- 5*setting21[['rct.size']]
setting21[['title']] <- 'Setting 21'

setting22 <- setting2
setting22[['rct.size']] <- 2*setting21[['rct.size']]
setting22[['title']] <- 'Setting 22'

setting23 <- setting2
setting23[['n.cov']] <- 10
setting23[['title']] <- 'Setting 23'

setting3 <- setting2
setting3[['title']] <- 'Setting 3'
setting3[['mixture.frac']] <- 1

setting31 <- setting3
setting31[['rwd.size']] <- 5*setting31[['rct.size']]
setting31[['title']] <- 'Setting 31'

setting32 <- setting3
setting32[['rct.size']] <- 2*setting32[['rct.size']]
setting32[['title']] <- 'Setting 32'

setting33 <- setting3
setting33[['n.cov']] <- 10
setting33[['title']] <- 'Setting 33'

setting4 <- setting2
setting4[['b.coef']]  <-  c(0, 0, 0, 0,  1)
setting4[['title']] <- 'Setting 4'

setting5 <- setting4
setting5[['b.coef']]  <-  c(0, 0, 0, 0.5, 0.5)
setting5[['title']] <- 'Setting 5'

setting6 <- setting4
setting6[['b.coef']]  <-  c(0, 0, 1/3, 1/3, 1/3)
setting6[['title']] <- 'Setting 6'

setting7 <- setting4
setting7[['b.coef']]  <-  c(0, 0.25, 0.25, 0.25, 0.25)
setting7[['title']] <- 'Setting 7'

runSimulation(setting=setting1, mainFolder = mainFolder)

runSimulation(setting=setting3, mainFolder = mainFolder)

