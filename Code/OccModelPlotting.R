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
pacman::p_load(pacman, tidyverse, here, jagsUI)

### Finish loading packages ###

### Load in Workspace from occ model run ###
load(workspace)

#############################################
############Model visualization##############
#############################################

#Assign model output to an object

mod <- "model_output"

### Trace plots ###
par(mfrow = c(), ask = F, mar = c())

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

