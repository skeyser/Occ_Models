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
jdf <- read.csv(file = here::here("Data_BBS/Generated DFs/jdf.csv"))
jdf <- jdf[, -1]
yrs <- read.csv(file = here::here("Data_BBs/Generated DFs/Years.Occ.csv"))
#Checkout the z.prime matrix 
dim(z.prime)

#z.prime is spp x sites x years
#Create empty matrix for handling 2-D matrices from z.prime
z2d.prime <- array2df(z.prime, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.5 <- array2df(z.prime.5, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.65 <- array2df(z.prime.65, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")
z2d.prime.75 <- array2df(z.prime.75, levels = list(spp.id = T, site.id = T, year.id = T), label.x = "Occupancy")

#UDF if col[x] is a factor make it an integer, if not make it a numeric
no.factors <- function(x){
  if(is.factor(x) == T){
    as.integer(x)
  } else {
    as.numeric(x)
  }
}

#Apply the UDF (no.factors)to all the DFs
z2d.prime[] <- lapply(z2d.prime, no.factors)
z2d.prime.5[] <- lapply(z2d.prime.5, no.factors)
z2d.prime.65[] <- lapply(z2d.prime.65, no.factors)
z2d.prime.75[] <- lapply(z2d.prime.75, no.factors)



#Combine z.dfs with info from spp.occ and jdf
df.prime <- z2d.prime %>% left_join( yrs, by = "year.id") %>%
            left_join(jdf, by = "site.id") %>% 
            left_join(spp.occ, by = "spp.id")

df.prime.5 <- z2d.prime.5 %>% left_join( yrs, by = "year.id") %>%
              left_join(jdf, by = "site.id") %>% 
              left_join(spp.occ, by = "spp.id")

df.prime.65 <- z2d.prime.65 %>% left_join( yrs, by = "year.id") %>%
               left_join(jdf, by = "site.id") %>% 
               left_join(spp.occ, by = "spp.id")

df.prime.75 <- z2d.prime.75 %>% left_join( yrs, by = "year.id") %>%
               left_join(jdf, by = "site.id") %>% 
               left_join(spp.occ, by = "spp.id")
