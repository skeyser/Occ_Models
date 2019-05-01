################################################
########Cleaning Data for Occupancy Models######
##########Author: Spencer R Keyser##############
################################################

#####Package Loading#####
library(pacman)
pacman::p_load(tidyverse, reshape2, ggplot2)

#####Packages Loaded######

#####Adding in UDFs for Occ######

expit <- function(p){
  return(exp(p) / (exp(p) + 1))
}

logit <- function(p, e = 0.0005){
  p = ifelse(p == 1, p- e, p)
  p = ifelse(p == 0, p + e, p)
  return(log(p / (1 - p)))
}

######Standardizing functions######
######Generic Function Gelman's Suggestion######

standardise <- function(xmat, stdev = 1, marg = c(1, 2, 3)){
  mean.xmat = mean(as.vector(xmat), na.rm = T)
  sd.xmat = sd(as.vector(xmat), na.rm = T)
  std.xmat = apply(xmat, marg, function(x){
    (x - mean.xmat) / (stdev * sd.xmat)})
  retun(std.xmat)
}

#####Standardizing Vectors#####
standardise.vector <- function(x, stdev = 1){(x - mean(x, na.rm = T))/
    (stdev * sd(x, na.rm = T))}

#####Ending Standardization Functions#####

#####Import species specific data#####