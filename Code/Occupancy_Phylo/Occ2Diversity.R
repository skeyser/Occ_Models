#########################################################################################
#######################Script for turning MSOM Output back into##########################
#######################Community Data for Biodiversity Analyses##########################
#########################################################################################
##########################Script Created By: Spencer R Keyser############################
##################################### 7/16/2019 #########################################
#########################################################################################

#Clean the environment
rm(list = ls())

#Load in Packages

library("pacman")
pacman::p_load("tidyverse", "plyr", "reshape2", "here", "arrayhelpers")

#Finish Package Loading 

#Load in Work spaces of z-matrices 
load(file = here::here("R Workspace/Output/Output4Analysis/Muscicapida.RData"))

#Checkout the z.prime matrix 
dim(z.prime)

#z.prime is spp x sites x years
#Create empty matrix for handling 2-D matrices from z.prime
z2d.prime <- array2df(z.prime, levels = list(Spp.id = T, Site.id = T, Year.id = T), label.x = "Occupancy")
z2d.prime.5 <- array2df(z.prime.5, levels = list(Spp.id = T, Site.id = T, Year.id = T), label.x = "Occupancy")
z2d.prime.65 <- array2df(z.prime.65, levels = list(Spp.id = T, Site.id = T, Year.id = T), label.x = "Occupancy")
z2d.prime.75 <- array2df(z.prime.75, levels = list(Spp.id = T, Site.id = T, Year.id = T), label.x = "Occupancy")
