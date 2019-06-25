######################################
####Script made for visualization#####
####And plotting of results from #####
####            MSOM             #####
#### Created by: Spencer R Keyser#####
########## 06/07/2019 ################
######################################

#Clean Workspace
rm(list = ls())

#Package Loading 
library(pacman)
pacman::p_load(pacman, tidyverse, here, jagsUI, MCMCvis)

### Finish loading packages ###

### Load in Workspace from occ model run ###
output <- readRDS("E:/MSI_Research/MCMC Output/PelicaniformesSummaryMat.rds")

#############################################
############Model visualization##############
#############################################

#Assign model output to an object

mod <- output

### Trace plots ###
par(mfrow = c(3, 3), ask = F, mar = c(3, 4, 2, 2))

traceplot(mod, parameters = c( 'int.p' ))
traceplot(mod, parameters = c( 'alpha' ))
traceplot(mod, parameters = c( 'p' ))
traceplot(mod, parameters = c( 'z' ))
traceplot(mod, parameters = c( 'int.psi' ))
traceplot(mod, parameters = c( 'eps.psi' ))
traceplot(mod, parameters = c( 'w.bcr' ))
traceplot(mod, parameters = c( 'epsID.psi' ))
traceplot(mod, parameters = c( 'sigma.psi' ))

###Finish trace plots###

